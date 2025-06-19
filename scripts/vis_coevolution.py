#!/usr/bin/env python3
"""
vis_coevolution.py

Visualize PLMC co-evolution couplings on a reference sequence.
Numbering is relative to the first sequence in the MSA; gaps are removed.
Produces both a scatter plot and a heatmap of coupling strengths.
"""
import argparse
import os
import sys
import csv
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np


def load_reference_mapping(msa_file):
    """
    Load the MSA and return a mapping from alignment columns (1-based)
    to reference positions (1-based, no gaps) for the first sequence.
    Also returns the reference length (ungapped).
    """
    aln = AlignIO.read(msa_file, "fasta")
    ref = aln[0].seq
    mapping = {}
    pos = 0
    for idx, aa in enumerate(ref, start=1):
        if aa != "-":
            pos += 1
            mapping[idx] = pos
    return mapping, pos


def load_couplings(csv_file):
    """
    Read a CSV of i,j,strength (1-based alignment positions);
    returns list of (i, j, strength) tuples.
    """
    couplings = []
    with open(csv_file) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            i = int(row['i'])
            j = int(row['j'])
            strength = float(row['strength'])
            couplings.append((i, j, strength))
    return couplings


def map_couplings(couplings, mapping):
    """
    Map alignment positions to reference positions; drop any pair
    where either alignment position is a gap in reference.
    """
    mapped = []
    for i, j, strength in couplings:
        if i in mapping and j in mapping:
            mapped.append((mapping[i], mapping[j], strength))
    return mapped


def plot_scatter(mapped, ref_length, out_scatter):
    """
    Create and save a scatter plot of couplings.
    """
    x = np.array([i for i, _, _ in mapped])
    y = np.array([j for _, j, _ in mapped])
    s = np.array([strength for _, _, strength in mapped])

    plt.figure(figsize=(8, 8))
    pts = plt.scatter(x, y, c=s, cmap='viridis', marker='o', s=20)
    plt.colorbar(pts, label='Coupling strength')
    plt.xlabel('Position (ref) i')
    plt.ylabel('Position (ref) j')
    plt.title('Co-evolution couplings (scatter)')
    plt.xlim(1, ref_length)
    plt.ylim(1, ref_length)
    plt.gca().set_aspect('equal', 'box')
    plt.tight_layout()
    plt.savefig(out_scatter)
    plt.close()


def plot_heatmap(mapped, ref_length, out_heatmap):
    """
    Create and save a heatmap of coupling strengths.
    """
    # Initialize matrix
    mat = np.zeros((ref_length, ref_length))
    for i, j, strength in mapped:
        mat[i-1, j-1] = strength
        mat[j-1, i-1] = strength  # symmetrical

    plt.figure(figsize=(8, 6))
    im = plt.imshow(mat, origin='lower', aspect='auto', interpolation='nearest')
    plt.colorbar(im, label='Coupling strength')
    plt.xlabel('Position (ref)')
    plt.ylabel('Position (ref)')
    plt.title('Co-evolution couplings (heatmap)')
    plt.tight_layout()
    plt.savefig(out_heatmap)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize PLMC co-evolution couplings')
    parser.add_argument('--msa',      required=True, help='MSA file (FASTA or .aln)')
    parser.add_argument('--couplings',required=True, help='CSV of i,j,strength')
    parser.add_argument('--out',      required=True, help='Output prefix for image files')
    args = parser.parse_args()

    if not os.path.exists(args.msa):
        sys.exit(f"MSA file not found: {args.msa}")
    if not os.path.exists(args.couplings):
        sys.exit(f"Couplings CSV not found: {args.couplings}")

    mapping, ref_len = load_reference_mapping(args.msa)
    couplings = load_couplings(args.couplings)
    mapped = map_couplings(couplings, mapping)

    if not mapped:
        sys.exit("No couplings map to non-gap reference positions; nothing to plot.")

    # Scatter plot
    scatter_file = f"{args.out}_scatter.png"
    plot_scatter(mapped, ref_len, scatter_file)
    print(f"Saved scatter plot to {scatter_file}")

    # Heatmap
    heatmap_file = f"{args.out}_heatmap.png"
    plot_heatmap(mapped, ref_len, heatmap_file)
    print(f"Saved heatmap to {heatmap_file}")

if __name__ == '__main__':
    main()
