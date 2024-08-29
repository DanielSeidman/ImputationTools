zcat ${1}.vcf.gz | awk 'BEGIN{
        OFS="\t";
}
($1~"^#"){
        print
}
($1!~"^#"){
        for(x=10;x<=NF;x++){
                split($x,genotype,":");
                if(genotype[1]=="."){
                        sub(/^.*:/,"./.:",$x)
                }
        };
        print
}' | bcftools view -Oz -o ${1}_for_beagle.vcf.gz
