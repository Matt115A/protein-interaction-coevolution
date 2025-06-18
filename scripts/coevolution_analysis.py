# scripts/coevolution_analysis.py
#!/usr/bin/env python3
"""
Concatenate the two MSAs and run plmc to compute co‐evolutionary couplings.
Produces a CSV of (i,j,strength).
"""
import argparse
import os
import subprocess
from Bio import AlignIO

def concat_alignments(dirpath):
    a = AlignIO.read(os.path.join(dirpath, "proteinA.aln"), "fasta")
    b = AlignIO.read(os.path.join(dirpath, "proteinB.aln"), "fasta")
    if len(a) != len(b):
        raise ValueError("Number of sequences differs between A and B")
    out_path = os.path.join(dirpath, "concat.fasta")
    with open(out_path, "w") as fh:
        for ra, rb in zip(a, b):
            fh.write(f">{ra.id}\n{str(ra.seq)+str(rb.seq)}\n")
    return out_path

def run_plmc(seqfile, params_out, couplings_out, plmc_cmd="plmc"):
    cmd = [
        plmc_cmd,
        "-o", params_out,
        "-c", couplings_out,
        "-le", "16",    # L2 regularization on fields
        "-l2", "0.01",  # L2 regularization on couplings
        "-seq", seqfile
    ]
    subprocess.run(cmd, check=True)

def parse_couplings(cpl_file, csv_file):
    with open(cpl_file) as fh, open(csv_file,"w") as out:
        out.write("i,j,strength\n")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            i,j,score = line.split()[:3]
            out.write(f"{i},{j},{score}\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--indir", required=True,
                   help="dir with proteinA.aln & proteinB.aln")
    p.add_argument("--out", required=True,
                   help="CSV of couplings")
    p.add_argument("--plmc-cmd", default="plmc",
                   help="plmc executable on your PATH")
    args = p.parse_args()

    params = args.out.replace(".csv", ".params")
    couplings = args.out.replace(".csv", ".couplings")

    seqfile = concat_alignments(args.indir)
    run_plmc(seqfile, params, couplings, args.plmc_cmd)
    parse_couplings(couplings, args.out)

if __name__ == "__main__":
    main()
