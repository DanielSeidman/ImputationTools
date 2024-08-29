module load bwa
module load perl
module load R

scaffoldName=$1
wgsvcf=$2
outputsummary=$3



#creates a sam for the reference of interest
bwa mem -t 2 $REFERENCE FSJbeadchipSeq.fq > FSJbeadchipSeq.sam

awk '(!($1 ~ "^@")){print $1,$2}' FSJbeadchipSeq.sam | sort > wgs.pair.txt

perl seqLocation.pl FSJbeadchipSeq.sam FSJbeadchipSeqLoc.txt

Rscript BeadchipProbes_mod.R finalFSJbeadchip.csv FSJbeadchipSeqLoc.txt FSJbeadchipSeqLocFSJgenomeV2_15May2023Fixed.txt


#creates a comparison sam for the beadchip file
bwa mem -t 2 final.assembly.fasta FSJbeadchipSeq.fq > comparison.sam

awk '(!($1 ~ "^@")){print $1,$2}' comparison.sam | sort > comparison.pair.txt

perl seqLocation.pl comparison.sam comparisonSeqLoc.txt

Rscript BeadchipProbes_mod.R finalFSJbeadchip.csv comparisonSeqLoc.txt comparison_15May2023Fixed.txt


#To handle inconsistent alleles between the two references, we use the following script
join <(sort wgs.pair.txt) <(sort comparison.pair.txt) | awk '($2!=$3){print}' | awk '{split($1,scaff,"caffold_"); print scaff[1]scaff[2]}' > mismatchedAlleleSNPsBetweenReferences.txt

join <(sort mismatchedAlleleSNPsBetweenReferences.txt) <(awk '{split($1,scaff,"caffold_"); print scaff[1]scaff[2], $2"p"$7}' FSJbeadchipSeqLocFSJgenomeV2_15May2023Fixed.txt | sort) | awk '{print $2}' > mismatchedAlleleSNPsBetweenReferences.remapped.txt

awk 'NR==FNR{snps[$1]=1}NR!=FNR{if($2 in snps){print FNR}}' mismatchedAlleleSNPsBetweenReferences.txt FSJfullPedFiltDogFINAL12July2016final.map > snpIndicesOfMismatches.txt

awk '
BEGIN{
	complements[0]=0; 
	complements["A"]="T"; 
	complements["T"]="A"; 
	complements["C"]="G"; 
	complements["G"]="C"
}
NR==FNR{
	indices[$1]=1
}
NR!=FNR{
	count=1;
	for(x=7;x<=NF;x+=2){
		if((x-5)/2 in indices){
			$x=complements[$x];
			$(x+1)=complements[$(x+1)]
		}
	};
 	print
}' snpIndicesOfMismatches.txt FSJfullPedFiltDogFINAL12July2016final.ped > FSJfullPedFiltDogFINAL12July2016final.partial_complemented.ped 

#awk 'NR==FNR{translation[$2]=$4}NR!=FNR{if($2 in translation){$1=translation[$2];$2=translation[$2]; print}}' idUpdate.txt FSJfullPedFiltDogFINAL12July2016final.partial_complemented.ped > FSJfullPedFiltDogFINAL12July2016final.rename.repositioned.ped
awk 'NR==FNR{translation[$2]=$4}NR!=FNR{if($2 in translation){$1=translation[$2];$2=translation[$2]; print}}' idUpdate.txt FSJfullPedFiltDogFINAL12July2016final.ped > FSJfullPedFiltDogFINAL12July2016final.rename.repositioned.ped


awk 'NR==FNR{split($1,scaff,"caffold_"); reposition[scaff[1]scaff[2]]=$7}NR!=FNR{$4=reposition[$2]; print}' FSJbeadchipSeqLocFSJgenomeV2_15May2023Fixed.txt FSJfullPedFiltDogFINAL12July2016final.map > FSJfullPedFiltDogFINAL12July2016final.rename.repositioned.map



#Plink is weird. Some versions don't work properly with some types of commands, hence the constant loading/unloading of versions. There is likely a better way to have done this.
module unload plink;module load plink/1.9
#cat FSJfullPedFiltDogFINAL12July2016finalNewPos.ped | sed 's/G/Cx/' | sed 's/C/Gx/' | sed 's/T/Ax/' | sed 's/A/Tx/' | sed 's/Cx/C/' | sed 's/Gx/G/' | sed 's/Ax/A/' | sed 's/Tx/T/' > FSJfullPedFiltDogFINAL12July2016finalNewPos2.ped
ln -s FSJfullPedFiltDogFINAL12July2016finalNewPos.map FSJfullPedFiltDogFINAL12July2016finalNewPos2.map

#chromsInGroup=$(zcat /scratch/nchen11_lab/dseidmanProcesses/collective_analysis/all_samples_${scaffoldName}.vcf.gz | awk '(! ($1 ~ /^#/)){print $1}' | sort | uniq | awk 'BEGIN{val=""}{if(val!=""){val=val","}val=val$1}END{print val}')

chromsInGroup=$scaffoldName

plink --file FSJfullPedFiltDogFINAL12July2016final.rename.repositioned --chr $chromsInGroup --make-bed --out beadchip-${scaffoldName} --chr-set 41 no-xy no-mt --output-chr chrMT

module unload plink; module load plink

plink2 --bfile beadchip-${scaffoldName} --set-all-var-ids @p# --recode vcf --out beadchipAltRef-${scaffoldName} --chr-set 41 no-xy no-mt --output-chr chrMT

awk '(!($1 ~ /^#/)){print $3}' beadchipAltRef-${scaffoldName}.vcf > beadChipVarsUpdate.txt

module unload plink; module load plink/1.9

#/scratch/nchen11_lab/dseidmanProcesses/collective_analysis/imputationDir/out.gt.vcf.gz
#../collective_analysis/all_samples_chrom_test_chr24.vcf.gz


plink --vcf $2 --chr $chromsInGroup --chr-set 41 no-xy no-mt --out evaluateMarkers-${scaffoldName} --export vcf --output-chr chrMT

module unload plink; module load plink

plink2 --vcf evaluateMarkers-${scaffoldName}.vcf --set-all-var-ids @p#  --chr-set 41 no-xy no-mt --out evaluateMarkers-${scaffoldName}-2 --export vcf --output-chr chrMT

plink2 --vcf evaluateMarkers-${scaffoldName}-2.vcf --chr-set 41 no-xy no-mt --extract beadChipVarsUpdate.txt --out evaluateMarkers-${scaffoldName}-3 --export vcf --output-chr chrMT

awk '(!($1 ~ /^#/)){print $3}' evaluateMarkers-${scaffoldName}-3.vcf > evalVarsUpdate.txt

plink2 --vcf beadchipAltRef-${scaffoldName}.vcf --chr-set 41 no-xy no-mt --extract evalVarsUpdate.txt --out beadchipAltRef-${scaffoldName}-2 --export vcf --output-chr chrMT


bash pythonOnlyConcordanceJoint.sh $scaffoldName evaluateMarkers-${scaffoldName}-2.vcf $scaffoldName


#module load gatk
#module load python
#module load bcftools

#awk 'NR==FNR{++a[$1]} !a[$1]' problemSampleH.txt <(bcftools query -l evaluateMarkers-${scaffoldName}-3.vcf) > vcfsamples.txt
#grep "H" vcfsamples.txt > hsamples.txt
#grep "L" vcfsamples.txt > lsamples.txt
#bcftools view -S ^problemSampleH.txt  beadchipAltRef-${scaffoldName}.vcf | bcftools view -S hsamples.txt - > beadchipAltRef-${scaffoldName}-sorted.vcf
#bcftools view -S ^problemSampleH.txt  evaluateMarkers-${scaffoldName}-3.vcf | bcftools view -S hsamples.txt - > evaluateMarkers-${scaffoldName}-3-sorted.vcf
#python /scratch/dseidman/concordanceChecker.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > hsummary_$3.txt

#bcftools view -S lsamples.txt beadchipAltRef-${scaffoldName}.vcf > beadchipAltRef-${scaffoldName}-sorted.vcf
#bcftools view -S lsamples.txt evaluateMarkers-${scaffoldName}-3.vcf > evaluateMarkers-${scaffoldName}-3-sorted.vcf
#python /scratch/dseidman/concordanceChecker.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > lsummary_$3.txt

#bcftools view -S ^problemSampleH.txt beadchipAltRef-${scaffoldName}.vcf | bcftools view -S vcfsamples.txt - > beadchipAltRef-${scaffoldName}-sorted.vcf
#bcftools view -S ^problemSampleH.txt evaluateMarkers-${scaffoldName}-3.vcf | bcftools view -S vcfsamples.txt - > evaluateMarkers-${scaffoldName}-3-sorted.vcf
#python /scratch/dseidman/concordanceChecker.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > summary_$3.txt

#gatk Concordance -R $REFERENCE --truth fsjWG-3-2-23Samples3-${scaffoldName}.vcf --eval FSJV3VCF4-${scaffoldName}.vcf --summary summary.tsv
