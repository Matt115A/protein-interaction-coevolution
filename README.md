# Protein Interaction Coevolution Pipeline

> A modular pipeline for detecting co-evolving residues between interacting protein pairs, inspired by [Hopf et al., 2021, Nature Communications](https://www.nature.com/articles/s41467-021-21636-z).

---

## Overview

This pipeline enables the detection and visualization of co-evolving positions between two protein sequences, leveraging large-scale homology search, multiple sequence alignment (MSA), and statistical coupling analysis. It is designed for flexibility, reproducibility, and ease of use.

---

## Features
- Automated BLAST homology search with assembly-level filtering for taxonomic diversity
- Pairing of homologs from the same genome assembly
- Concatenated MSA construction and quality filtering
- Coevolution analysis using PLMC
- Publication-quality heatmaps and network visualizations
- Modular, script-driven workflow

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone <repo-url>
   cd protein-interaction-coevolution
   ```
2. **Run setup:**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```
   This will install all required Python dependencies automatically.
3. **Install MAFFT:**
   - MAFFT is required for multiple sequence alignment. Install via your package manager (e.g., `brew install mafft` on macOS, `apt-get install mafft` on Ubuntu) or from [the official website](https://mafft.cbrc.jp/alignment/software/).

4. **Install PLMC:**
   - PLMC is required for coevolution analysis. Download and compile PLMC from [the official repository](https://github.com/debbiemarkslab/plmc):
     ```bash
     git clone https://github.com/debbiemarkslab/plmc.git
     cd plmc
     make
     ```
   - After compilation, copy the `plmc` binary to a directory in your `PATH` (e.g., `/usr/local/bin`) **or** place it in the root of this repository. The pipeline will look for `plmc` in your `PATH` or in the current working directory.
   - To test your installation:
     ```bash
     plmc -h
     ```

---

## Quick Start

1. **Prepare your input:**
   - Place your two query FASTA files as `queryA.fasta` and `queryB.fasta` in `inputs/find_homologues/PROJECT/`.
   - Set the `project` name in `project.yaml`.

2. **Run the pipeline:**
   ```bash
   make find        # Find and pair homologs
   make msa         # Build and filter concatenated MSA
   make coevolution # Run PLMC coevolution analysis
   make heatmap     # Generate heatmaps and network visualizations
   ```

---

## Pipeline Steps

### 1. `make find`
- **Description:**
  - Runs `scripts/blast_and_pair.py` to perform BLASTP searches for both query proteins against NCBI nr, filters by assembly, and pairs homologs from the same genome.
- **Input:**
  - `inputs/find_homologues/PROJECT/queryA.fasta`
  - `inputs/find_homologues/PROJECT/queryB.fasta`
- **Output:**
  - `results/find_homologues/PROJECT/paired_sequences.json`

### 2. `make msa`
- **Description:**
  - Runs `scripts/build_msas.sh` to construct a concatenated MSA of paired homologs using MAFFT, with post-processing to remove low-quality columns.
- **Input:**
  - `results/find_homologues/PROJECT/paired_sequences.json`
- **Output:**
  - `results/msas/PROJECT/proteinAB.aln`

### 3. `make coevolution`
- **Description:**
  - Runs `scripts/coevolution_analysis.py` to compute coevolutionary couplings using PLMC, outputting a CSV of residue-residue coupling strengths.
- **Input:**
  - `results/msas/PROJECT/proteinAB.aln`
- **Output:**
  - `results/coevolution/PROJECT/coevolution_results.csv`

### 4. `make heatmap`
- **Description:**
  - Runs `scripts/generate_heatmaps.py` to create heatmaps and network diagrams of the strongest coevolving residue pairs.
- **Input:**
  - All previous outputs
- **Output:**
  - `results/heatmaps/PROJECT/` (PNG images, CSV mappings)

---

## Output
- **Paired homologs:** JSON file with all paired sequences
- **MSA:** Concatenated and filtered alignment in FASTA format
- **Coevolution results:** CSV of coupling scores
- **Visualizations:**
  - Global heatmaps at multiple thresholds
  - Network diagrams of inter-protein couplings
  - Mapped coupling CSVs for downstream analysis

---

## Dependencies
- Python 3.7+
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [scipy](https://scipy.org/)
- [pyyaml](https://pyyaml.org/)
- [tqdm](https://tqdm.github.io/)
- [pysam](https://pysam.readthedocs.io/en/latest/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (external, for MSA)
- [PLMC](https://github.com/debbiemarkslab/plmc) (external, for coevolution analysis)

---

## Reference
This pipeline is inspired by:

Hopf, T. A., Ingraham, J. B., Poelwijk, F. J., Schärfe, C. P. I., Springer, M., Sander, C., & Marks, D. S. (2021). Mutation effects predicted from sequence co-variation. *Nature Communications*, 12, 1232. [https://www.nature.com/articles/s41467-021-21636-z](https://www.nature.com/articles/s41467-021-21636-z)

If you use this pipeline, please cite the above work and this repository.

---

## License

Distributed under the MIT License. See `LICENSE` for details.
