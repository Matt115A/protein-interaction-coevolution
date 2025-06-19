#!/usr/bin/env bash
# scripts/build_msas.sh
# Usage: build_msas.sh paired_sequences.json output_directory
set -euo pipefail

PAIRED_JSON="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

python3 - "$PAIRED_JSON" "$OUTDIR" <<'EOF'
import sys, json, random

paired_json = sys.argv[1]
outdir      = sys.argv[2]
pairs       = json.load(open(paired_json))

if not pairs:
    sys.exit("ERROR: no pairs in JSON")

# baseline = the first entry, always included
first   = pairs[0]
lenA0   = len(first["proteinA"]["seq"])
lenB0   = len(first["proteinB"]["seq"])
minA    = int(0.8 * lenA0)
maxA    = int(1.2 * lenA0)
minB    = int(0.8 * lenB0)
maxB    = int(1.2 * lenB0)

# collect one eligible pair per assembly (excluding the baseline asm)
by_asm = {}
for p in pairs[1:]:
    asm = p["assembly"]
    if asm == first["assembly"] or asm in by_asm:
        continue
    seqA = p["proteinA"]["seq"]
    seqB = p["proteinB"]["seq"]
    if not (minA <= len(seqA) <= maxA and minB <= len(seqB) <= maxB):
        continue
    by_asm[asm] = p

# shuffle candidates
candidates = list(by_asm.values())
random.shuffle(candidates)

# build final list, ensuring no identical concat-seq is picked twice
picked      = [first]
seen_concat = { first["proteinA"]["seq"] + first["proteinB"]["seq"] }

for p in candidates:
    if len(picked) >= 100:
        break
    combo = p["proteinA"]["seq"] + p["proteinB"]["seq"]
    if combo in seen_concat:
        continue
    seen_concat.add(combo)
    picked.append(p)

# write FASTA
with open(f"{outdir}/proteinAB.fasta","w") as fa:
    for p in picked:
        asm  = p["assembly"]
        idA  = p["proteinA"]["id"]
        idB  = p["proteinB"]["id"]
        seqA = p["proteinA"]["seq"]
        seqB = p["proteinB"]["seq"]
        header = f">{asm}|{idA}|{idB}"
        fa.write(f"{header}\n{seqA}{seqB}\n")
EOF

# Build MSA with MAFFT on the concatenated sequences
mafft --auto "$OUTDIR"/proteinAB.fasta > "$OUTDIR"/proteinAB.aln.tmp

# Post‐process: drop columns gapped in the ref AND <80% non‐gapped overall
python3 - "$OUTDIR"/proteinAB.aln.tmp "$OUTDIR"/proteinAB.aln <<'EOF'
import sys
from Bio import AlignIO

in_aln, out_aln = sys.argv[1], sys.argv[2]
aln = AlignIO.read(in_aln, "fasta")
num_seqs = len(aln)
threshold = 0.8
ref_seq = aln[0].seq  # first record

# find columns to keep
keep_cols = []
for col_idx in range(aln.get_alignment_length()):
    if ref_seq[col_idx] == "-":
        non_gaps = sum(rec.seq[col_idx] != "-" for rec in aln)
        if non_gaps / num_seqs < threshold:
            continue
    keep_cols.append(col_idx)

# write filtered alignment
with open(out_aln, "w") as out:
    for rec in aln:
        filtered = "".join(rec.seq[i] for i in keep_cols)
        out.write(f">{rec.id}\n{filtered}\n")
EOF

# clean up temporary
rm "$OUTDIR"/proteinAB.aln.tmp
