#!/bin/bash
#SBATCH --partition=rosalind  --time 5-00:00:00 --mem=300g -c 12
#SBATCH --mail-user=dseidman@ur.rochester.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd /scratch/nchen11_lab/dseidmanProcesses/collective_analyses

module load java



java -Djava.io.tmpdir=/scratch/dseidman/tempdir -jar /scratch/dseidman/bin/beagle.22Jul22.46e.jar gt=../beadchipAnalysis/beadchip-out.vcf.gz out=beadchip_phased burnin=3 iterations=12
