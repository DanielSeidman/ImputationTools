rm -f markerFile24.txt

join <(awk '{for(x=1;x<=NF;x++){print $x}}' /scratch/nchen11_lab/linkageMap/chr24/chr24.CRIMAP2507.frameworkSNPs | sort) <(awk '{print $1,$2}' /scratch/nchen11_lab/linkageMap/chr24/chr24.loc | sort) | awk '{print $1}' > /scratch/nchen11_lab/dseidmanProcesses/collective_analysis/imputationDir/gigi24frameworkSNPs.txt

awk 'NR==FNR{ind1[$1]=$2} NR!=FNR{print ind1[$1]}' /scratch/nchen11_lab/linkageMap/chr24/chr24.loc gigi24frameworkSNPs.txt > gigi24frameworkSNPs.toID.txt

awk '  function abs(v) {return v < 0 ? -v : v} BEGIN{lastval=0}NR==FNR{split($1,snpname,"caffold_"); snps[snpname[1]snpname[2]]=$3}NR!=FNR{if($1 in snps){if(lastval!=0){print abs(snps[$1]-lastval)}; lastval=snps[$1]}}' HifiSal_Dec22_beadchipSeqLoc_mappositions.txt gigi24frameworkSNPs.toID.txt | awk 'BEGIN{strToPrint="map marker dist "}{strToPrint=strToPrint" "$1}END{print strToPrint}' >> markerFile24.txt

printf "\n" >> markerFile24.txt

awk '{print $2,"chr"$1"p"$4}' /scratch/nchen11_lab/dseidmanProcesses/beadchipAnalysis/FSJfullPedFiltDogFINAL12July2016final.rename.repositioned.map > oldSnpToTempName24.txt

join <(sort gigi24frameworkSNPs.toID.txt) <(sort oldSnpToTempName24.txt) | awk '{print $2}' > snpNamesForFreqs24.txt

bcftools view ../all_samples_chrom_test_chr24.vcf.gz -R targets24.txt > ../all_samples_chrom_test_chr24_scaffoldSNPsonly.vcf

plink2 --vcf ../all_samples_chrom_test_chr24_scaffoldSNPsonly.vcf --chr-set 41 no-xy no-mt --output-chr chrMT --set-all-var-ids @p# --freq --out chr24

awk '($1 !~ "^#"){print "set markers "(NR-1)"\tallele freqs "(1-$5)" "$5}' chr24.afreq >> markerFile24.txt

printf "\n" >> markerFile24.txt

printf "set markers 51 data\n\n" >> markerFile24.txt


awk '($1=="#CHROM"){for(x=10;x<=NF;x++){strings[x]=$x};fieldNum=NF}($1 !~ "^#"){for(x=10;x<=NF;x++){split($x,qual,":|/|\\|"); strings[x]=strings[x]" "(qual[1]+1)" "(qual[2]+1)}}END{for(x=10;x<=NF;x++){if(strings[x] !~ ","){print strings[x]}}}' ../all_samples_chrom_test_chr24_scaffoldSNPsonly.vcf >> markerFile24.txt
