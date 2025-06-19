#!/usr/bin/env python3
"""
scripts/coevolution_analysis.py

Concatenate MSAs if needed, then run PLMC to compute co‐evolutionary couplings.
Produces a CSV of (i,j,strength).
Supports either two separate ALN files or a single concatenated MSA (proteinAB.aln).
"""
import argparse
import os
import sys
import logging
from Bio import AlignIO
from subprocess import Popen, PIPE, STDOUT

def concat_alignments(dirpath):
    """
    Concatenate proteinA.aln + proteinB.aln into concat.fasta
    """
    a_path = os.path.join(dirpath, "proteinA.aln")
    b_path = os.path.join(dirpath, "proteinB.aln")
    logging.info(f"Concatenating {a_path} + {b_path}")
    a = AlignIO.read(a_path, "fasta")
    b = AlignIO.read(b_path, "fasta")
    if len(a) != len(b):
        raise ValueError("Number of sequences differs between A and B")
    out_path = os.path.join(dirpath, "concat.fasta")
    with open(out_path, "w") as fh:
        for ra, rb in zip(a, b):
            fh.write(f">{ra.id}\n{ra.seq}{rb.seq}\n")
    return out_path

def get_sequence_file(indir):
    """
    Return the path to the MSA file: prefer single proteinAB.aln, else concat two alignments.
    """
    ab_path = os.path.join(indir, "proteinAB.aln")
    if os.path.exists(ab_path):
        logging.info(f"Using existing concatenated MSA: {ab_path}")
        return ab_path
    return concat_alignments(indir)

def run_plmc(seqfile, params_out, couplings_out, plmc_cmd="plmc", maxiter=5):
    """
    Run PLMC on the given sequence file with a maximum number of iterations.
    Streams PLMC stdout/stderr live to your logger.
    Returns:
      0 if PLMC ran to completion,
      1 if it hit max-iterations (non-fatal),
    Exits on any other error.
    """
    cmd = [
        plmc_cmd,
        "-o", params_out,
        "-c", couplings_out,
        "-le", "16",
        "-l2", "0.01",
        "-t", "0.01",
        "-m", str(maxiter),
        seqfile
    ]
    logging.info(f"Running: {' '.join(cmd)}")
    proc = Popen(cmd, stdout=PIPE, stderr=STDOUT, text=True)

    # Stream PLMC output line by line
    for line in proc.stdout:
        logging.info(line.rstrip())

    proc.wait()
    ret = proc.returncode

    if ret == 0:
        return 0

    # PLMC uses exit code 1 when hitting max iterations but still writes output
    if ret == 1 and "MAXIMUMITERATION" in open(params_out, errors='ignore').read():
        logging.warning("PLMC hit maximum iterations (-m %d); continuing with whatever it produced", maxiter)
        return 1

    if ret < 0:
        sig = -ret
        sys.exit(f"PLMC terminated by signal {sig} (e.g. SIGSEGV). "
                 "This often means too few effective sequences in your MSA.")

    # Any other non-zero is a failure
    sys.exit(f"PLMC failed (exit code {ret}). See log above for details.")

def parse_couplings(cpl_file, csv_file):
    """
    Parse PLMC couplings file into CSV of i,j,strength.
    """
    with open(cpl_file) as fh, open(csv_file, "w") as out:
        out.write("i,j,strength\n")
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            i, j, score = line.split()[:3]
            out.write(f"{i},{j},{score}\n")

def main():
    parser = argparse.ArgumentParser(description='PLMC co-evolution analysis')
    parser.add_argument("--indir",  required=True,
                        help="Directory with MSA: proteinAB.aln or proteinA.aln+proteinB.aln")
    parser.add_argument("--out",    required=True,
                        help="Output CSV of couplings")
    parser.add_argument("--plmc-cmd", default="plmc",
                        help="plmc executable")
    parser.add_argument("-m", "--maxiter", type=int, default=50,
                        help="Maximum L-BFGS iterations for PLMC")
    args = parser.parse_args()

    # Ensure output directory exists
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Determine input MSA
    seqfile = get_sequence_file(args.indir)

    # Prepare output filenames
    base, _ = os.path.splitext(args.out)
    params_file    = base + ".params"
    couplings_file = base + ".couplings"

    # Run PLMC
    status = run_plmc(seqfile,
                      params_out=params_file,
                      couplings_out=couplings_file,
                      plmc_cmd=args.plmc_cmd,
                      maxiter=args.maxiter)

    # Validate output
    if not os.path.exists(couplings_file) or os.path.getsize(couplings_file) < 10:
        if status == 1:
            sys.exit(
                f"PLMC hit the maximum iterations (-m {args.maxiter}) but did not produce a valid couplings file.\n"
                "Try increasing --maxiter, using --fast, or improving your MSA diversity."
            )
        else:
            sys.exit(f"PLMC output {couplings_file} is missing or too small")

    # Parse results
    parse_couplings(couplings_file, args.out)
    logging.info(f"Wrote couplings CSV to {args.out}")

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr
    )
    main()
