To use this repo:

chmod +x setup.sh
./setup.sh

Overall, this code allows you to test two sequences for co-evolving positions and surfaces. It is still evolving code!


# protien-interaction-coevolution

commands: 

	make find
	make msa
	make coevolution
	make heatmap

Instructions: 

## find

Place your two inputs in inputs/find_homologues/PROJECT/ queryA.fasta and queryB.fasta. Define PROJECT in the project.yaml.

	make find

## msa

	make msa

## coevolution

Runs plmc with default parameters

	make coevolution

## heatmap

Generates heatmaps at different connection thresholds. Also generates between-protein connection maps at different thresholds. 

	make heatmap
