#!/bin/bash
### Environment variables 
export TMP=${WORKDIR}/output/tmp ;
export REF=${WORKDIR}/ref/Nipponbare_chr.fasta ;
export INPUT=${WORKDIR}/output/BAM ;
export BAM=${WORKDIR}/output/BAM;
export VCF=${WORKDIR}/output/VCF;

mkdir -p $BAM $VCF;

export CORE=1 ;
val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1)) ;
fi

## Read one Sample Per Rank 
cp $WORKDIR/scripts/Phase1/Phase1.txt $WORKDIR/scripts/Phase2/Phase2.prefix.txt ;
DATA=`sed -n ${LINE}p  $WORKDIR/scripts/Phase2/Phase2.prefix.txt`;
PREFIX=`basename $DATA .1.fastq.gz` ;

echo "--------- Step 1: sam-sort----------------------------------------------------------"
echo "samtools sort -T /demo/tmp/${PREFIX} -o $BAM/${PREFIX}.sorted.bam $BAM/${PREFIX}.sorted.fxmt.mkdup.addrep.bam" ;
samtools sort -T $TMP/${PREFIX} -o $BAM/${PREFIX}.sorted.bam $BAM/${PREFIX}.sorted.fxmt.mkdup.addrep.bam ;
echo "--------- Step 2: sam-index --------------------------------------------------------"
echo "samtools index $BAM/${PREFIX}.sorted.bam" ;
samtools index $BAM/${PREFIX}.sorted.bam ;
echo "--------- Step 3: HaplotypeCaller ---------------------------------------------------"
echo "gatk --java-options "-Xmx2g -Xms2g" HaplotypeCaller --input $BAM/${PREFIX}.sorted.bam --output $VCF/${PREFIX}.snps.indels.g.vcf.gz --reference $REF --emit-ref-confidence GVCF --min-base-quality-score 20 --bam-output $BAM/${PREFIX}.assembledHC.bam --tmp-dir /demo/tmp/${PREFIX}" ;
mkdir -p $TMP/${PREFIX} ;
gatk --java-options "-Xmx2g -Xms2g" HaplotypeCaller --input $BAM/${PREFIX}.sorted.bam --output $VCF/${PREFIX}.snps.indels.g.vcf.gz --reference $REF --emit-ref-confidence GVCF --min-base-quality-score 20 --bam-output $BAM/${PREFIX}.assembledHC.bam --tmp-dir $TMP/${PREFIX} ;
echo "---------- Step 4: Tabix ---------------------------------------------------------------"  
echo "tabix -f $VCF/$PREFIX.snps.indels.g.vcf.gz" ;
tabix -f $VCF/$PREFIX.snps.indels.g.vcf.gz ;
