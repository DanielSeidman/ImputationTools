scaffoldName=$1
wgsvcf=$2
outputsummary=$3

module load gatk
module load python
module load bcftools

awk '
BEGIN{
	OFS="\t";
	FS="\t";
	t["A"]="T";
	t["T"]="A";
	t["C"]="G";
	t["G"]="C"
}
NR==FNR{
	refs[$3]=$4;
	alts[$3]=$5
}
(NR!=FNR && $1 !~ "^#"){
	status="anomalous";
	if(refs[$3]==$4 && (alts[$3]=="." || alts[$3]==$5)){
		status="correct"
	}else if(((alts[$3]=="." && refs[$3]!=t[$5]) || alts[$3]==$4) && refs[$3]==$5){
		status="ref-alt-flip"
	}else if(refs[$3]==t[$4] && (alts[$3]=="." || alts[$3]==t[$5])){
		status="reverse-comp"
	}else if((alts[$3]==t[$4] || alts[$3]==".") && refs[$3]==t[$5]){
		status="rev-and-flip"
	};
	print $3,refs[$3],alts[$3],$4,$5,status
}' evaluateMarkers-${scaffoldName}-3.vcf  beadchipAltRef-${scaffoldName}.vcf > flipStatus.txt

awk '
BEGIN{
	OFS="\t";
	t["A"]="T";
        t["T"]="A";
        t["C"]="G";
        t["G"]="C"
}
NR==FNR{
	flipState[$1]=$6; 
}
(NR!=FNR && (($1 ~ "^#") || (flipState[$3]=="correct"))){
	print
}
(NR!=FNR && ($1 !~ "^#") && (flipState[$3]!="correct")){
	if(flipState[$3]=="rev-and-flip" || flipState[$3]=="ref-alt-flip"){
		a=$4
		$4=$5;
		$5=a;
	};
	if(flipState[$3]=="rev-and-flip" || flipState[$3]=="reverse-comp"){
		$4=t[$4];
		$5=t[$5];
	}
        flipSNP["0/0"]="1/1";
	flipSNP["1/1"]="0/0";
	flipSNP["./."]="./.";
	flipSNP["0/1"]="0/1";
	flipSNP["1/0"]="1/0";
	if(flipState[$3]=="ref-alt-flip" || flipState[$3]=="rev-and-flip"){
		for(x=10;x<=NF;x++){
			$x=flipSNP[$x];
	}
	};
 	print
}' flipStatus.txt beadchipAltRef-${scaffoldName}-2.vcf > beadchipAltRef2-${scaffoldName}.vcf

awk 'NR==FNR{++a[$1]} !a[$1]' problemSampleH.txt <(bcftools query -l evaluateMarkers-${scaffoldName}-3.vcf) > vcfsamples.txt
grep "H" vcfsamples.txt > hsamples.txt
grep "L" vcfsamples.txt > lsamples.txt
bcftools view -S ^problemSampleH.txt  beadchipAltRef2-${scaffoldName}.vcf | bcftools view -S hsamples.txt - > beadchipAltRef-${scaffoldName}-sorted.vcf
bcftools view -S ^problemSampleH.txt  evaluateMarkers-${scaffoldName}-3.vcf | bcftools view -S hsamples.txt - > evaluateMarkers-${scaffoldName}-3-sorted.vcf
python /scratch/dseidman/concordanceChecker2.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > hsummary_$3.txt

bcftools view -S lsamples.txt beadchipAltRef2-${scaffoldName}.vcf > beadchipAltRef-${scaffoldName}-sorted.vcf
bcftools view -S lsamples.txt evaluateMarkers-${scaffoldName}-3.vcf > evaluateMarkers-${scaffoldName}-3-sorted.vcf
python /scratch/dseidman/concordanceChecker2.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > lsummary_$3.txt

bcftools view -S ^problemSampleH.txt beadchipAltRef2-${scaffoldName}.vcf | bcftools view -S vcfsamples.txt - > beadchipAltRef-${scaffoldName}-sorted.vcf
bcftools view -S ^problemSampleH.txt evaluateMarkers-${scaffoldName}-3.vcf | bcftools view -S vcfsamples.txt - > evaluateMarkers-${scaffoldName}-3-sorted.vcf
python /scratch/dseidman/concordanceChecker2.py evaluateMarkers-${scaffoldName}-3-sorted.vcf beadchipAltRef-${scaffoldName}-sorted.vcf > summary_$3.txt
