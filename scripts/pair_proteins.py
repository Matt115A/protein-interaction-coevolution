# scripts/pair_proteins.py
#!/usr/bin/env python3
"""
Given a raw multi‐FASTA, group sequences by assembly and extract two target proteins per genome.
Outputs JSON like:
[ { "assembly": "...",
    "proteinA": {"id":"...","seq":"..."},
    "proteinB": {"id":"...","seq":"..."} },
  … ]
"""
import argparse, json
from collections import defaultdict
from Bio import SeqIO

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in", dest="in_fasta", required=True)
    p.add_argument("--protA-name", required=True,
                   help="substring to identify protein A (e.g. 'GeneX')")
    p.add_argument("--protB-name", required=True,
                   help="substring to identify protein B (e.g. 'GeneY')")
    p.add_argument("--out", required=True, help="paired JSON")
    args = p.parse_args()

    by_asm = defaultdict(list)
    for rec in SeqIO.parse(args.in_fasta, "fasta"):
        # Expect header like >assembly|protein_id
        asm, pid = rec.id.split("|",1)
        by_asm[asm].append((pid, str(rec.seq)))

    pairs = []
    for asm, seqs in by_asm.items():
        a = [s for s in seqs if args.protA_name in s[0]]
        b = [s for s in seqs if args.protB_name in s[0]]
        if len(a)==1 and len(b)==1:
            pairs.append({
                "assembly": asm,
                "proteinA": {"id": a[0][0], "seq": a[0][1]},
                "proteinB": {"id": b[0][0], "seq": b[0][1]},
            })

    with open(args.out, "w") as fh:
        json.dump(pairs, fh, indent=2)

if __name__ == "__main__":
    main()
