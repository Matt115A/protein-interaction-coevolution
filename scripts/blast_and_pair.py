#!/usr/bin/env python3
"""
scripts/blast_and_pair.py

1. BLAST two query proteins via NCBIWWW.qblast (remote), restricted to mammals.
2. Map each hit accession → its assemblies via Entrez.elink.
3. For each assembly, form all combinations of A- and B-hits (cartesian product).
4. Fetch the full-length protein sequences in batches.
5. Write paired_sequences.json.

Optimizations:
- Uses qblast (no external process spawning).
- Parses XML in-memory via StringIO.
- Captures only first HSP per alignment up to n_hits.
- DEBUG logs for command and BLAST runtime.
"""
import argparse
import json
import os
import sys
import time
from io import StringIO
from collections import defaultdict
from itertools import product
import logging
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW

# Configure logging
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr
    )


def blast_hits(sequence, n_hits, tax_filter):
    logging.info(f"Starting qblast for top {n_hits} hits, filter={tax_filter}")
    start = time.time()
    result_handle = NCBIWWW.qblast(
        program='blastp',
        database='nr',
        sequence=sequence,
        expect=0.00001,
        hitlist_size=n_hits,
        entrez_query=tax_filter
    )
    xml_data = result_handle.read()
    result_handle.close()
    duration = time.time() - start
    logging.info(f"qblast finished in {duration:.1f}s, parsing XML data")

    blast_record = NCBIXML.read(StringIO(xml_data))
    accs = []
    seen = set()
    for alignment in blast_record.alignments:
        if len(accs) >= n_hits:
            break
        acc = alignment.accession
        if acc in seen:
            continue
        # only capture first HSP
        accs.append(acc)
        seen.add(acc)
    logging.info(f"Collected {len(accs)} accessions from BLAST")
    return accs


def map_to_assemblies(accession, email):
    Entrez.email = email
    try:
        handle = Entrez.elink(
            dbfrom='protein', db='assembly', id=accession,
            linkname='protein_assembly'
        )
        recs = Entrez.read(handle)
        handle.close()
        asms = [link['Id']
                for ls in recs
                for db in ls.get('LinkSetDb', [])
                for link in db.get('Link', [])]
        logging.debug(f"{accession}→assemblies: {asms}")
        return asms
    except Exception as e:
        logging.warning(f"Error mapping {accession}: {e}")
        return []


def fetch_sequences(accessions, email, batch_size=100):
    logging.info(f"Fetching {len(accessions)} sequences in batches of {batch_size}")
    Entrez.email = email
    seqs = {}
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        logging.debug(f"Fetching batch {i}-{i+len(batch)-1}")
        handle = Entrez.efetch(
            db='protein', id=','.join(batch),
            rettype='fasta', retmode='text'
        )
        for rec in SeqIO.parse(handle, 'fasta'):
            seqs[rec.id] = rec
        handle.close()
        time.sleep(0.4)
    logging.info(f"Fetched total {len(seqs)} sequences")
    return seqs


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description='Find and pair homologues by genome via qblast+Entrez'
    )
    parser.add_argument('--queryA', required=True, help='protein A FASTA')
    parser.add_argument('--queryB', required=True, help='protein B FASTA')
    parser.add_argument('--n-hits', type=int, default=50, help='top N BLAST hits')
    parser.add_argument('--email', required=True, help='Entrez email')
    parser.add_argument('--out', required=True, help='output JSON path')
    parser.add_argument('--tax-filter', default='mammalia[Organism]', help='Entrez taxonomy filter')
    args = parser.parse_args()

    logging.info('--- Pipeline start ---')
    # Load query sequences
    seqA = SeqIO.read(args.queryA, 'fasta')
    seqB = SeqIO.read(args.queryB, 'fasta')

    # BLAST queries
    hitsA = blast_hits(str(seqA.seq), args.n_hits, args.tax_filter)
    hitsB = blast_hits(str(seqB.seq), args.n_hits, args.tax_filter)

    # Map hits to assemblies
    asm_map = defaultdict(lambda: {'A': [], 'B': []})
    for acc in hitsA:
        for asm in map_to_assemblies(acc, args.email):
            asm_map[asm]['A'].append(acc)
    for acc in hitsB:
        for asm in map_to_assemblies(acc, args.email):
            asm_map[asm]['B'].append(acc)

    # Filter assemblies with at least one of each
    assemblies = [a for a,m in asm_map.items() if m['A'] and m['B']]
    logging.info(f"Assemblies with both hits: {len(assemblies)}")

    # Fetch unique sequences
    all_accs = sorted({*hitsA, *hitsB})
    seq_records = fetch_sequences(all_accs, args.email)

    # Build pairs
    pairs = []
    for asm in assemblies:
        for a_acc, b_acc in product(asm_map[asm]['A'], asm_map[asm]['B']):
            recA = seq_records.get(a_acc)
            recB = seq_records.get(b_acc)
            if recA and recB:
                pairs.append({
                    'assembly': asm,
                    'proteinA': {'id': recA.id, 'seq': str(recA.seq)},
                    'proteinB': {'id': recB.id, 'seq': str(recB.seq)}
                })
    # Write output
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as out_f:
        json.dump(pairs, out_f, indent=2)
    logging.info(f"Wrote {len(pairs)} pairs to {args.out}")
    logging.info('--- Pipeline end ---')

if __name__ == '__main__':
    main()
