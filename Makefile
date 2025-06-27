PROJECT := protein-interaction-coevolution

find:
	python scripts/blast_and_pair.py \
	  --queryA inputs/find_homologues/$(PROJECT)/queryA.fasta \
	  --queryB inputs/find_homologues/$(PROJECT)/queryB.fasta \
	  --n-hits 500 \
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

heatmap:
	python3 scripts/generate_heatmaps.py --project $(PROJECT)
