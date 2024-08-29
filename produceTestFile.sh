activename=${1}
module load vcftools

zcat ${activename}.vcf.gz | awk -v "gqThresh=$2" '
BEGIN{
        OFS="\t"
}
($1~"^#"){
        print
}
($1=="#CHROM"){
        for(x=10;x<=NF;x++){
                samples[x]=$x
        }
}
($1!~"^#"){
        for(x=10;x<=NF;x++){
                split($x,genotype,":");
		if(genotype[1]=="."){
			sub(/^.*:/,"./.:",$x)
		}
                if((genotype[4]<gqThresh)){
                        sub(/^.*:/,"./.:",$x)
                };
        };
        print
}' | gzip > ${activename}_${2}.vcf.gz

activename=${activename}_${2}

module load vcftools

vcftools --gzvcf ${activename}.vcf.gz --max-missing $4 --stdout --recode | gzip > ${activename}_${4}.vcf.gz



activename=${activename}_${4}

comparisonname=$activename

zcat ${activename}.vcf.gz | awk 'BEGIN{
        OFS="\t";
}
($1~"^#"){
        print
}
($1=="#CHROM"){
        for(x=10;x<=NF;x++){
                samples[x]=$x
        }
}
($1!~"^#"){
        for(x=10;x<=NF;x++){
                split($x,genotype,":");
                if(genotype[1]=="."){
                        sub(/^.*:/,"./.:",$x)
                }
        };
        print
}' | gzip > ${activename}_for_beagle.vcf.gz


zcat ${activename}.vcf.gz | awk -v "downFrac=$3" 'BEGIN{
        OFS="\t";
}
($1~"^#"){
        print
}
($1=="#CHROM"){
        for(x=10;x<=NF;x++){
                samples[x]=$x
        }
}
($1!~"^#"){
        for(x=10;x<=NF;x++){
                split($x,genotype,":");
		if(genotype[1]=="."){
                        sub(/^.*:/,"./.:",$x)
                }
                if(((samples[x]~"^H") && x%2==0 && rand()>(1-downFrac))){
                        sub(/^.*:/,"./.:",$x)
                };
        };
        print
}' | gzip > ${activename}_downsample${3}.vcf.gz
activename=${activename}_downsample${3}

bash imputationDir/alphaImputeTests/convertToAlphaImputeGenotypesSpecific.sh ${comparisonname}.vcf.gz ${2}_${4}

bash imputationDir/alphaImputeTests/convertToAlphaImputeGenotypesSpecific.sh ${activename}.vcf.gz $2_$3

awk -v "modVal=$5" '{toPrint=$1; for(x=2;x<=NF;x++){if(x%modVal==0){ toPrint=toPrint" "$x}}; print toPrint}' Genotypes_specific_$2_${3}.txt > Genotypes_specific_$2_${3}_mod${5}.txt
awk -v "modVal=$5" '{toPrint=$1; for(x=2;x<=NF;x++){if(x%modVal==0){ toPrint=toPrint" "$x}}; print toPrint}' Genotypes_specific_$2_${4}.txt > Genotypes_specific_$2_${4}_mod${5}.txt
cp Genotypes.txt imputationDir/alphaImputeTests/.

cd imputationDir/alphaImputeTests/

sbatch alphaImputeRun.sh
