To use this repo:

chmod +x setup.sh
./setup.sh

# protien-interaction-coevolution

commands: 

	make find
	make msa
	make coevolution

Instructions: 

## find

Place your two inputs in inputs/find_homologues/PROJECT/ queryA.fasta and queryB.fasta. Define PROJECT in the project.yaml.

	make find

## msa

	make msa

## coevolution

Runs plmc with default parameters

	make coevolution
