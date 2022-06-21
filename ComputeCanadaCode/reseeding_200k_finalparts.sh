#!/bin/bash

#SBATCH --job-name=reseeding_finalparts 
#SBATCH --output=outputgraph/output_reseeding_finalparts.txt 
#SBATCH --error=outputgraph/error_reseeding_finalparts.txt
#SBATCH --mail-user=ariel.kiara@gmail.com
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=REQUEUE 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=186G
#SBATCH --time=0-65:00

module load arch/avx2 StdEnv/2016.4

module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

Rscript Reseeding200k_finalparts_computecanadaversion.R 
