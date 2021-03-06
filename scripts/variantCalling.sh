#! /bin/bash

## MToolBox calling ##

CONF="./variantCalling.conf"

MToolBox.sh -i $CONF > log.txt

### Varscan2 uses realigned mtDNA bam file from MToolBox to call the variants 

BAM_FILE=$1 #Bam list file
SAMPLE_ID=${BAM_FILE##*/}
SAMPLE_ID=${SAMPLE_ID%.bam}

VARSCAN_OUT_DIR="./iPSC_Analysis/Varscan_OUTPUT/BaseQual30"

if [ ! -d $VARSCAN_OUT_DIR ]
then
        mkdir $VARSCAN_OUT_DIR
else
        echo "$VARSCAN_OUT_DIR exists"   
fi

SCRIPT_PATH="./Scripts"
TARGETS="mito.bed"
REF_FILE_MT="./rCRS.fa"
Q="30"
echo "Samtools starts to make pileup file, then Varscan starts to call the variants"
#Coverage#
COV_DIR="${VARSCAN_OUT_DIR}/Coverage"

if [ ! -d $COV_DIR ]
then
        mkdir $COV_DIR
else
        echo "$COV_DIR exists"   
fi
##Coverge per base on mitochrondrial geonome##
PILEUP_rCRS="$COV_DIR/${SAMPLE_ID}.pileup"
samtools mpileup -d 0 -q $Q -Q $Q -l $TARGETS -f $REF_FILE_MT $BAM_FILE > $PILEUP_rCRS

###Varscan calling#######

echo 'Starting Varscan2'

VCF_FILErCRS="${VARSCAN_OUT_DIR}/${SAMPLE_ID}_snps.vcf"
java -Xmx12g -jar ${VSP}` mpileup2snp $PILEUP_rCRS --min-coverage 1 --min-reads2 1 --min-avg-qual $Q --min-var-freq 0.001  --output-vcf > $VCF_FILErCRS

echo 'VarScan variant calling has finished'

echo 'rm pileup file'
rm $PILEUP_rCRS

echo "Done."

## haplogroup defining using haploGrep2 ##
## only homoplasmies used to define mitochrondrial haplogroups ##

INPUT_VCF='Hip12765.haplogroup.vcf'
haploGrep2="./haplogrep-2.1.1.jar"

java -jar ${haploGrep2}  --format vcf --in $INPUT_VCF --out ${INPUT_VCF}.haplogrep --phylotree 17

