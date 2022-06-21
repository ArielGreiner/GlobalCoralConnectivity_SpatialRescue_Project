#!/bin/bash

#SBATCH --job-name=largestcluster_graphingtry1 
#SBATCH --output=outputgraph/output_largestcluster_graphingtry1.txt 
#SBATCH --error=outputgraph/error_largestcluster_graphingtry1.txt
#SBATCH --mail-user=ariel.kiara@gmail.com
#SBATCH --mail-type=BEGIN 
#SBATCH --mail-type=END 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=REQUEUE 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=96G
#SBATCH --time=0-16:00

module load arch/avx2 StdEnv/2016.4

module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0 

Rscript Plotting200k_partials/Plotting200k_largestcluster.R 
