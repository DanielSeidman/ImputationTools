#!/bin/bash
#SBATCH --partition=standard  --time 5-00:00:00 --mem=50g -c 8
#SBATCH --mail-user=dseidman@ur.rochester.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
REFERENCE=/scratch/nchen11_lab/FSJgenome/HifiSal_MAP_Dec22.fasta

cd /scratch/nchen11_lab/dseidmanProcesses/collective_analysis/imputationDir/alphaImputeTests


./AlphaImpute_Linux

