#!/bin/bash
## Variables
export TMP=${WORKDIR}/output/tmp ;
export REF=${WORKDIR}/ref/Nipponbare_chr.fasta ;
export gVCF=${WORKDIR}/output/gVCF;
export SNPs=${WORKDIR}/output/SNPs;
export INDELs=${WORKDIR}/output/INDELs;

mkdir -p $TMP $gVCF $SNPs $INDELs ;

## Rank based value 

val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1)) ;
fi

## Read one Split value Per Rank 
ChrName=`sed -n ${LINE}p $WORKDIR/scripts/Phase3/split.txt | awk '{print $1}'` ;
size=`sed -n ${LINE}p $WORKDIR/scripts/Phase3/split.txt | awk '{print $2}'` ;
Start=`sed -n ${LINE}p $WORKDIR/scripts/Phase3/split.txt | awk '{print $3}'` ;
End=`sed -n ${LINE}p $WORKDIR/scripts/Phase3/split.txt | awk '{print $4}'` ;
INPUT=`cat $WORKDIR/scripts/Phase3/Phase3.input.list`

## Step 1
echo " -----Step1, GATK: CombineGVCFs --------------"
gatk --java-options "-Xmx2g -Xms2g" CombineGVCFs --variant $INPUT --reference $REF --intervals $ChrName:$Start-$End --output ${gVCF}/Combine.$ChrName.$size.vcf.gz --tmp-dir $TMP ;
## Step 2
echo " ----- Step 2, GATK: GenotypeGVCFs -----------"
gatk --java-options "-Xmx2g -Xms2g" GenotypeGVCFs --variant $gVCF/Combine.$ChrName.$size.vcf.gz --reference $REF --intervals $ChrName:$Start-$End --output $gVCF/Genotype.$ChrName.$size.vcf.gz --tmp-dir $TMP ;
## Step 3
echo " ----- Step 3, GATK: SelectVariants (SNPs) -----------"
gatk --java-options "-Xmx2g -Xms2g" SelectVariants --variant $gVCF/Genotype.$ChrName.$size.vcf.gz --reference $REF -select-type SNP --output $SNPs/$ChrName.$size.vcf.gz --tmp-dir $TMP ;
## Step 4
echo " ----- Step 4, GATK: VariantFiltration (SNPs) --------"
gatk --java-options "-Xmx2g -Xms2g" VariantFiltration --variant $SNPs/$ChrName.$size.vcf.gz --reference $REF --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 20.0 || MQRankSum < -3.0 || ReadPosRanKSum < -3.0 || DP < 5.0" --filter-name snp_filter --output $SNPs/filtered_snps.$ChrName.$size.vcf.gz --tmp-dir $TMP ;
## Step 5
echo " ----- Step 5, GATK: SelectVariants (INDELs) -----------"
gatk --java-options "-Xmx2g -Xms2g" SelectVariants --variant $gVCF/Genotype.$ChrName.$size.vcf.gz --reference $REF -select-type INDEL --output $INDELs/$ChrName.$size.vcf.gz --tmp-dir $TMP;
## Step 6:
echo " ----- Step 6, GATK: VariantFiltration (INDELs) --------"
gatk --java-options "-Xmx2g -Xms2g" VariantFiltration --variant $INDELs/$ChrName.$size.vcf.gz --reference $REF --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 20.0 || MQRankSum < -3.0 || ReadPosRanKSum < -3.0 || DP < 5.0" --filter-name indel_filter --output $INDELs/filtered_indels.$ChrName.$size.vcf.gz --tmp-dir $TMP;
