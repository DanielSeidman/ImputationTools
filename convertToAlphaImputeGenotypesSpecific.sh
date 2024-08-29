
zcat $1 | awk 'BEGIN{conversion["./."]=3;conversion[".|."]=3;conversion["."]=3;conversion["0/0"]=0;conversion["0/1"]=1;conversion["1/1"]=2;conversion["0|1"]=1;conversion["1|0"]=1;conversion["0|0"]=0;conversion["1|1"]=2;}($1!~"^#"){strToPrint=$1"p"$2;for(x=10;x<=NF;x++){split($x,genotype,":");g=genotype[1]; strToPrint=strToPrint" "conversion[g]}; print strToPrint}' > genotypeMidpoint1.txt

module load bcftools

bcftools query -l $1 | awk '{split($1,sample,"_"); print sample[1]}'> vcfsamples.txt

awk 'NR==FNR{genotypes[FNR]=$1}NR!=FNR{for(x=2;x<=NF;x++){genotypes[(x-1)]=genotypes[(x-1)]" "$x}}END{for(x=1;x<=length(genotypes);x++){print samples[x]" "genotypes[x]}}' vcfsamples.txt genotypeMidpoint1.txt > Genotypes_specific_${2}.txt

#awk '($1~"^H"){print}($1~"^L"){for(x=2;x<=NF;x++){if($x==1){$x=3}};print}' Genotypes_no_missing.txt > Genotypes.txt

cp Genotypes_specific_${2}.txt Genotypes.txt
