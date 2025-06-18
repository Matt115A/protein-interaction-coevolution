# scripts/fetch_proteins.py
#!/usr/bin/env python3
"""
Fetch all protein sequences linked to a list of NCBI Assembly accessions.
Writes a single multi‐FASTA of all proteins to stdout or to --out.
"""
import argparse
import time
from Bio import Entrez, SeqIO

def fetch_protein_ids(assembly_id):
    """Use ELink to find all protein UIDs for a given assembly."""
    handle = Entrez.elink(dbfrom="assembly", db="protein",
                          id=assembly_id, linkname="assembly_protein")
    record = Entrez.read(handle)
    handle.close()
    uids = []
    for linkset in record:
        for link in linkset.get("LinkSetDb", []):
            uids += [l["Id"] for l in link["Link"]]
    return uids

def batch_fetch_fasta(uids, email, out_handle, batch_size=200):
    Entrez.email = email
    for i in range(0, len(uids), batch_size):
        chunk = uids[i:i+batch_size]
        handle = Entrez.efetch(db="protein", id=",".join(chunk),
                               rettype="fasta", retmode="text")
        out_handle.write(handle.read())
        handle.close()
        time.sleep(0.4)  # be gentle on NCBI

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genome-list", required=True,
                   help="one assembly accession per line")
    p.add_argument("--email", required=True, help="NCBI Entrez email")
    p.add_argument("--out", default="-",
                   help="output FASTA (default stdout)")
    args = p.parse_args()

    if args.out == "-":
        out_handle = sys.stdout
    else:
        out_handle = open(args.out, "w")

    all_uids = []
    with open(args.genome_list) as fh:
        for line in fh:
            asm = line.strip()
            if not asm or asm.startswith("#"):
                continue
            uids = fetch_protein_ids(asm)
            print(f"Found {len(uids)} proteins for {asm}", file=sys.stderr)
            all_uids.extend(uids)

    batch_fetch_fasta(all_uids, args.email, out_handle)

    if out_handle is not sys.stdout:
        out_handle.close()

if __name__ == "__main__":
    import sys
    main()
