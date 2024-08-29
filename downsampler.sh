#zcat all_samples_chrom_test3_chr24_bqsr_snps.vcf.gz | awk '
#BEGIN{
#	OFS="\t"
#}
#($1=="#CHROM"){
#	for(x=1;x<=NF;x++){
#		samples[x]=$x
#	}
#}
#($1~"^#"){
#	print
#}
#($1!~"^#"){
#	for(x=10;x<=NF;x++){
#		split($x,genotype,":"); 
#		if(samples[x]~"^H" && rand()>0.65){
#			sub(genotype[1],"./.",$x)
#		}
#	};
# 	print
#}' | gzip > all_samples_chrom_test3_chr24_bqsr_snps_downsample35.vcf.gz


zcat ${1}.vcf.gz | awk -v "downFrac=$2" 'BEGIN{
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
		if(((samples[x]~"^H") && x%2==0 && rand()>(1-downFrac))){
			sub(/^.*:/,"./.:",$x)
		};
	};
	print	
}' | gzip > ${1}_downsampled${2}.vcf.gz


