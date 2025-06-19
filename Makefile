PROJECT := protein-interaction-coevolution

find:
	python scripts/blast_and_pair.py \
	  --queryA inputs/find_homologues/$(PROJECT)/queryA.fasta \
	  --queryB inputs/find_homologues/$(PROJECT)/queryB.fasta \
	  --n-hits 5000 \
	  --email mp957@cam.ac.uk \
	  --out results/find_homologues/$(PROJECT)/paired_sequences.json

msa:
	bash scripts/build_msas.sh \
	  results/find_homologues/$(PROJECT)/paired_sequences.json \
	  results/msas/$(PROJECT)/

coevolution:
	python scripts/coevolution_analysis.py \
	  --indir results/msas/$(PROJECT)/ \
	  --out results/coevolution/$(PROJECT)/coevolution_results.csv


visualize: results/vis_coev_scatter.png results/vis_coev_heatmap.png

results/vis_coev_scatter.png results/vis_coev_heatmap.png: \
    results/proteinAB.aln \
    results/coevolution_results.csv \
    vis_coevolution.py
	@echo " ⤷ Generating co-evolution plots…"
	python3 vis_coevolution.py \
	  --msa  $< \
	  --couplings results/coevolution_results.csv \
	  --out  results/vis_coev