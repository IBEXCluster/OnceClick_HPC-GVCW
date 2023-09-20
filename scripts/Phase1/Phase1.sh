#!/bin/bash
### Environment variables 
export REF=${WORKDIR}/ref/Nipponbare_chr.fasta ;
export INPUT=${WORKDIR}/input ;
export BAM=${WORKDIR}/output/BAM ;
export TMP=${WORKDIR}/output/tmp ;
export CORE=1 ;

mkdir -p $TMP ;

val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1))
fi

## Read one Sample Per Rank 
DATA=`sed -n ${LINE}p  $WORKDIR/scripts/Phase1/Phase1.txt`;
PREFIX=`basename $DATA .1.fastq.gz` ;
SAMPLE=${PREFIX%*_*};


## Random wait before creating $BAM directory (to safeguard /luster file system)
export MAXWAIT=5
sleep $[ ( $RANDOM % $MAXWAIT )  + 1 ]s
mkdir -p $BAM;


## Step 1
echo "------------- Step 1 executing: BWA MEM ---------------------------------------"
echo "bwa mem -M -t $CORE $REF $INPUT/${PREFIX}.1.fastq.gz $INPUT/${PREFIX}.2.fastq.gz | samtools view -@ $CORE -b -S -h -q 30 - | samtools sort -T $TMP/${PREFIX} - > $BAM/$PREFIX.sorted.bam"
bwa mem -M -t $CORE $REF $INPUT/${PREFIX}.1.fastq.gz $INPUT/${PREFIX}.2.fastq.gz | samtools view -@ $CORE -b -S -h -q 30 - | samtools sort -T $TMP/${PREFIX} - > $BAM/$PREFIX.sorted.bam 
## Step 2
echo "------------- Step 2 executing: FixMateInformation  ------------------"
echo "gatk --java-options -Xmx2g -Xms2g FixMateInformation --INPUT $BAM/$PREFIX.sorted.bam --SORT_ORDER coordinate --OUTPUT $BAM/$PREFIX.sorted.fxmt.bam --TMP_DIR $TMP/${PREFIX}"
gatk --java-options "-Xmx2g -Xms2g" FixMateInformation --INPUT $BAM/$PREFIX.sorted.bam --SORT_ORDER coordinate --OUTPUT $BAM/$PREFIX.sorted.fxmt.bam --TMP_DIR $TMP/${PREFIX}
## Step 3
echo "------------ Step 3 executing: MarkDuplicates ------------"
echo "gatk --java-options -Xmx2g -Xms2g MarkDuplicates --INPUT $BAM/$PREFIX.sorted.fxmt.bam --METRICS_FILE $BAM/$PREFIX.metrics --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --TMP_DIR $TMP/${PREFIX}"
gatk --java-options "-Xmx2g -Xms2g" MarkDuplicates --INPUT $BAM/$PREFIX.sorted.fxmt.bam --METRICS_FILE $BAM/$PREFIX.metrics --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --TMP_DIR $TMP/${PREFIX}
## Step 4
echo "------------ Step 4 executing: AddOrReplaceReadGroups -------------------"
echo "gatk --java-options -Xmx2g -Xms2g AddOrReplaceReadGroups --INPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.addrep.bam --RGID $SAMPLE --RGPL Illumina --RGSM $SAMPLE --RGLB $SAMPLE --RGPU unit1 --RGCN BGI --SORT_ORDER coordinate --TMP_DIR $TMP/${PREFIX}"
gatk --java-options "-Xmx2g -Xms2g" AddOrReplaceReadGroups --INPUT $BAM/$PREFIX.sorted.fxmt.mkdup.bam --OUTPUT $BAM/$PREFIX.sorted.fxmt.mkdup.addrep.bam --RGID $SAMPLE --RGPL Illumina --RGSM $SAMPLE --RGLB $SAMPLE --RGPU unit1 --RGCN BGI --SORT_ORDER coordinate --TMP_DIR $TMP/${PREFIX}
