# ImputationTools

Many of these tools assume you have the absolute location of your reference sequence stored in an environmental variable, $REFERENCE.
You can set this with:

REFERENCE="/example/absolute/path/to/reference.fasta"
You can have this automatically set every time you start a cluster session by modifying your .bashrc file to include this, but is not necessary.


File locations of key data input:
/scratch/nchen11_lab/dseidmanProcesses/beadchipAnalysis/final.assembly.fasta
Nancy will also have this somewhere, but this is the theoretical final assembly that matches the beadchip snps, and is necessary if starting a new mapping of beadchip snps to a new reference.

/scratch/nchen11_lab/dseidmanProcesses/beadchipAnalysis/FSJbeadchipSeq.fq
Nancy will also have this somewhere, but this is the fasta for the beadchip calls I used for analyses

/scratch/nchen11_lab/dseidmanProcesses/beadchipAnalysis/FSJfullPedFiltDogFINAL12July2016final.map
I do not know why this says "dog" in it. I'm assuming I got part of that from the naming convention of the file I copied this from


Tool1: concordanceJoint.sh
Script for comparing a new set of calls to beadchip, accounting for variable orientation of scaffolds, in terms of directionality (3'-5' or 5'-3') or reference-alternate allele swap.
Takes three arguments:
scaffoldname, the name of the scaffold you are interested in running this with. 
wgsvcf, a location for the vcf of interest to use. The name comes from early attempts I was using where I thought I needed to include invariant sites to match, but now it assumes only a variant vcf. The match difficulties were a side effect of originally not properly reorienting scaffolds vs beadchip snps, which has since been rectified in this version.
outputsummary, a suffix for your output

Note: some midpoint names are hardcoded.


Tool2: BeadchipProbes_mod.R
Similar to what Nancy uses for beadchip probe processing, but with a few formatting modifications to function in my scripts.



Instructions for using other tools:

Alphaimpute was, at the time, obtainable through:

wget https://alphagenes.roslin.ed.ac.uk/wp/?ddownload=447
unzip AlphaImpute.zip

And more recently,

git clone https://github.com/AlphaGenes/AlphaImpute2.git
unzip AlphaImpute2.zip -d AlphaImpute2
cd AlphaImpute2
python setup.py build



I have a version of the former at: /scratch/dseidman/bin/AlphaImpute
and the latter at: /scratch/dseidman/bin/AlphaImpute2

But it may be wiser to re-download or have CIRC help install it



