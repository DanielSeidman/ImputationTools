#!/bin/bash
#SBATCH --partition=standard  --time 5-00:00:00 --mem=300g -c 12
#SBATCH --mail-user=dseidman@ur.rochester.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
REFERENCE=/scratch/nchen11_lab/FSJgenome/HifiSal_MAP_Dec22.fasta

cd /scratch/nchen11_lab/dseidmanProcesses/collective_analysis/imputationDir/

module load java


java -Djava.io.tmpdir=/scratch/dseidman/tempdir -jar /scratch/dseidman/bin/beagle.r1399.jar gt=all_samples_chrom_test3_chr24_bqsr_snps_biallelic_2_half_downsampled0.33.vcf.gz out=all_samples_chrom_test3_chr24_bqsr_snps_biallelic_2_half_downsampled0.33_beagle.vcf.gz > out.txt 2>error.txt
