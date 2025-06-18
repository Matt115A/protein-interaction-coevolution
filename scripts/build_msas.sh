#!/usr/bin/env bash
# scripts/build_msas.sh
# Usage: build_msas.sh paired_sequences.json output_directory
set -euo pipefail

PAIRED_JSON="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

# Split into two FASTAs
python3 - <<EOF
import sys, json
pairs = json.load(open(sys.argv[1]))
faA = open(f"{sys.argv[2]}/proteinA.fasta","w")
faB = open(f"{sys.argv[2]}/proteinB.fasta","w")
for p in pairs:
    asm = p["assembly"]
    faA.write(f">{asm}|{p['proteinA']['id']}\n{p['proteinA']['seq']}\n")
    faB.write(f">{asm}|{p['proteinB']['id']}\n{p['proteinB']['seq']}\n")
faA.close()
faB.close()
EOF

# Build MSAs with MAFFT (you can swap for clustalo etc.)
mafft --auto "$OUTDIR"/proteinA.fasta > "$OUTDIR"/proteinA.aln
mafft --auto "$OUTDIR"/proteinB.fasta > "$OUTDIR"/proteinB.aln
