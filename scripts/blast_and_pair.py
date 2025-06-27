#!/usr/bin/env python3
"""
scripts/blast_and_pair.py

Full end-to-end script with checkpointing, detailed logging,
robust assembly mapping, and continuous feedback.
Outputs diagnostic files and paired_sequences.json.
"""
import argparse
import json
import os
import sys
import time
import logging
from io import StringIO
from collections import defaultdict

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW

# Ensure project root is on PYTHONPATH so we can import utils
HERE = os.path.abspath(os.path.dirname(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, PROJECT_ROOT)

# Import the debug helper
from utils.debug_entrez import debug_entrez

# Cache for accession→assemblies
CACHE_MAP = {}
# UID mapping log handle
UID_LOG_FH = None


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr
    )


def blast_hits(fasta_path, n_hits, tax_filter, hit_file):
    """
    Run or load BLAST hits for a query FASTA.
    Saves to hit_file or loads if exists.
    """
    if os.path.exists(hit_file):
        logging.info(f"Loading existing hits from {hit_file}")
        with open(hit_file) as f:
            hits = json.load(f)
        logging.info(f"Loaded {len(hits)} hits for {fasta_path}")
        return hits
    seq = SeqIO.read(fasta_path, 'fasta')
    logging.info(f"Running qblast for {fasta_path}")
    start = time.time()
    handle = NCBIWWW.qblast(
        program='blastp', database='nr', sequence=str(seq.seq),
        expect=1e-5, hitlist_size=n_hits, entrez_query=tax_filter
    )
    xml_data = handle.read()
    handle.close()
    logging.info(f"qblast finished in {time.time()-start:.1f}s")
    rec = NCBIXML.read(StringIO(xml_data))
    hits, seen = [], set()
    for aln in rec.alignments:
        if len(hits) >= n_hits:
            break
        acc = aln.accession
        if acc in seen:
            continue
        hsp = aln.hsps[0]
        hits.append({
            'accession': acc,
            'evalue': hsp.expect,
            'bit_score': hsp.bits,
            'align_len': hsp.align_length
        })
        seen.add(acc)
    with open(hit_file, 'w') as f:
        json.dump(hits, f, indent=2)
    logging.info(f"Saved {len(hits)} hits to {hit_file}")
    return hits


def map_accession(accession, email):
    """
    Map protein accession → assembly accessions with full debug and a
    secondary nuccore→assembly fallback.
    """
    global CACHE_MAP, UID_LOG_FH
    if accession in CACHE_MAP:
        logging.debug(f"CACHE HIT for {accession}: {CACHE_MAP[accession]}")
        return CACHE_MAP[accession]

    Entrez.email = email
    logging.info(f"Mapping accession {accession} to assemblies")
    start = time.time()
    assemblies = []

    # --- 1) ESEARCH: accession → protein UID ---
    uid = None
    try:
        es = Entrez.esearch(db='protein', term=f"{accession}[Accession]")
        rec = Entrez.read(es)
        es.close()
        logging.debug(f"esearch record for {accession}: {rec}")
        ids = rec.get('IdList', [])
        if ids:
            uid = ids[0]
            logging.info(f"  → Found protein UID: {uid}")
        else:
            logging.warning(f"  → No protein UID found for {accession}")
    except Exception as e:
        logging.error(f"  esearch ERROR for {accession}: {e} (skipping UID mapping)")

    UID_LOG_FH.write(f"{accession}\t{uid or ''}\n")
    UID_LOG_FH.flush()

    # --- 2) PRIMARY ELink: protein → assembly ---
    asm_uids = []
    if uid:
        try:
            el = Entrez.elink(dbfrom='protein', db='assembly', id=uid)
            recs = Entrez.read(el)
            el.close()
            logging.debug(f"elink protein→assembly raw: {recs}")
            for linkset in recs:
                for linkdb in linkset.get('LinkSetDb', []):
                    for link in linkdb.get('Link', []):
                        asm_uids.append(link['Id'])
            logging.info(f"  → PRIMARY elink returned {len(asm_uids)} assembly UIDs")
        except Exception as e:
            logging.error(f"  elink ERROR for UID {uid}: {e} (continuing fallback)")
    else:
        logging.info("  → Skipping primary elink (no UID)")

    # --- 3) SECONDARY FALLBACK: protein → nuccore → assembly ---
    if not asm_uids and uid:
        logging.info("  → No direct assemblies: trying protein→nuccore→assembly fallback")
        nuc_ids = []
        try:
            el2 = Entrez.elink(dbfrom='protein', db='nuccore', id=uid)
            rec2 = Entrez.read(el2)
            el2.close()
            for ls in rec2:
                for ldb in ls.get('LinkSetDb', []):
                    for link in ldb.get('Link', []):
                        nuc_ids.append(link['Id'])
            logging.info(f"    → Found {len(nuc_ids)} nuccore IDs")
        except Exception as e:
            logging.error(f"    protein→nuccore elink ERROR: {e}")

        if nuc_ids:
            try:
                el3 = Entrez.elink(dbfrom='nuccore', db='assembly', id=','.join(nuc_ids))
                rec3 = Entrez.read(el3)
                el3.close()
                for ls in rec3:
                    for ldb in ls.get('LinkSetDb', []):
                        for link in ldb.get('Link', []):
                            asm_uids.append(link['Id'])
                logging.info(f"    → Fallback elink returned {len(asm_uids)} assembly UIDs")
            except Exception as e:
                logging.error(f"    nuccore→assembly elink ERROR: {e}")
        else:
            logging.info("    → Skipping nuccore→assembly (no nuccore IDs)")

    # --- 4) ESUMMARY if we got any assembly UIDs ---
    if asm_uids:
        try:
            sum_h = Entrez.esummary(db='assembly', id=','.join(asm_uids))
            summary = Entrez.read(sum_h)
            sum_h.close()
            docs = summary.get('DocumentSummarySet', {}).get('DocumentSummary', [])
            logging.info(f"  → esummary returned {len(docs)} docs")
            for doc in docs:
                acc = doc.get('AssemblyAccession')
                if acc:
                    assemblies.append(acc)
            logging.info(f"  → Parsed {len(assemblies)} AssemblyAccession values")
        except Exception as e:
            logging.error(f"  esummary ERROR: {e}")
    else:
        logging.info("  → No assembly UIDs to summary")

    # --- 5) XML fallback (xref) as last resort ---
    if not assemblies and uid:
        try:
            xf = Entrez.efetch(db='protein', id=uid, retmode='xml')
            docs = Entrez.read(xf)
            xf.close()
            for doc in docs:
                xrefs = doc.get('GBSeq_xref', [])
                for xr in xrefs:
                    if xr.get('GBXref_dbname') == 'Assembly' and xr.get('GBXref_id'):
                        assemblies.append(xr['GBXref_id'])
            logging.info(f"  → XML fallback extracted {len(assemblies)} assemblies")
        except Exception as e:
            logging.error(f"  XML fallback ERROR for {accession}: {e}")

    logging.info(f"Mapped in {time.time()-start:.2f}s → {len(assemblies)} assemblies: {assemblies}")
    CACHE_MAP[accession] = assemblies
    return assemblies


def fetch_sequences(accessions, email, batch_size):
    Entrez.email = email
    seqs = {}
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        handle = Entrez.efetch(db='protein', id=','.join(batch), rettype='fasta', retmode='text')
        for rec in SeqIO.parse(handle, 'fasta'):
            base = rec.id.split('.')[0]
            seqs[base] = rec
        handle.close()
        time.sleep(0.4)
    return seqs


def main():
    setup_logging()
    parser = argparse.ArgumentParser(description='BLAST & pair homologues')
    parser.add_argument('--queryA', required=True)
    parser.add_argument('--queryB', required=True)
    parser.add_argument('--n-hits', type=int, default=50)
    parser.add_argument('--email', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--tax-filter', default='mammalia[Organism]')
    args = parser.parse_args()

    # Verify Entrez availability before proceeding
    debug_entrez(email=args.email)

    out_dir = os.path.dirname(args.out)
    os.makedirs(out_dir, exist_ok=True)
    hitsA_file = os.path.join(out_dir,'hitsA_list.json')
    hitsB_file = os.path.join(out_dir,'hitsB_list.json')
    uid_file    = os.path.join(out_dir,'uid_mapping.tsv')
    summary_file= os.path.join(out_dir,'assembly_summary.tsv')

    global UID_LOG_FH
    UID_LOG_FH = open(uid_file,'a') if os.path.exists(uid_file) else open(uid_file,'w')
    if os.path.getsize(uid_file)==0:
        UID_LOG_FH.write('accession\tuid\n')

    logging.info('--- Pipeline start ---')
    # BLAST hits
    hitsA = blast_hits(args.queryA, args.n_hits, args.tax_filter, hitsA_file)
    hitsB = blast_hits(args.queryB, args.n_hits, args.tax_filter, hitsB_file)

    # Map
    asm_map = defaultdict(lambda: {'A': [], 'B': []})
    logging.info(f"Mapping {len(hitsA)} QueryA hits...")
    for idx, hit in enumerate(hitsA, 1):
        acc = hit['accession']
        logging.info(f"[{idx}/{len(hitsA)}] {acc}")
        asms = map_accession(acc, args.email)
        for a in asms:
            asm_map[a]['A'].append(acc)
    logging.info(f"Mapping {len(hitsB)} QueryB hits...")
    for idx, hit in enumerate(hitsB, 1):
        acc = hit['accession']
        logging.info(f"[{idx}/{len(hitsB)}] {acc}")
        asms = map_accession(acc, args.email)
        for a in asms:
            asm_map[a]['B'].append(acc)

    # Summary
    with open(summary_file, 'w') as f:
        f.write('assembly\tA_hits\tB_hits\n')
        for a, counts in asm_map.items():
            f.write(f"{a}\t{len(counts['A'])}\t{len(counts['B'])}\n")
    logging.info(f"Wrote assembly summary to {summary_file}")
    assemblies = [a for a, c in asm_map.items() if c['A'] and c['B']]
    logging.info(f"Assemblies with both hits: {len(assemblies)}")

    # Fetch sequences
    all_accs = sorted({h['accession'] for h in hitsA + hitsB})
    seqs = fetch_sequences(all_accs, args.email, 100)

    # Pair
    pairs = []
    for a in assemblies:
        for x in asm_map[a]['A']:
            for y in asm_map[a]['B']:
                r1, r2 = seqs.get(x), seqs.get(y)
                if r1 and r2:
                    pairs.append({
                        'assembly': a,
                        'proteinA': {'id': r1.id, 'seq': str(r1.seq)},
                        'proteinB': {'id': r2.id, 'seq': str(r2.seq)}
                    })
    with open(args.out, 'w') as f:
        json.dump(pairs, f, indent=2)
    logging.info(f"Wrote {len(pairs)} pairs to {args.out}")

    UID_LOG_FH.close()
    logging.info('--- Pipeline end ---')


if __name__ == '__main__':
    main()
