#!/bin/bash
## Variables
export TMP=${WORKDIR}/output/tmp ;
export REF=${WORKDIR}/ref/Nipponbare_chr.fasta ;
export OUTPUT=${WORKDIR}/output ;

mkdir -p $OUTPUT/SNPs_per_Chr $OUTPUT/INDELs_per_Chr ;
mkdir -p $OUTPUT/SNPs_Var_Table $OUTPUT/INDELs_Var_Table ;

## Rank based value 
val=$1;
if [ $val -eq 0 ]
then 
  LINE=1;
else
 LINE=$((val + 1)) ;
fi

## Read one Split value Per Rank 
CHR=`sed -n ${LINE}p ${WORKDIR}/scripts/Phase4/distribution.txt | awk '{print $1}'` ;
CHUNK=`sed -n ${LINE}p ${WORKDIR}/scripts/Phase4/distribution.txt | awk '{print $2}'` ;

set INPUT_SNP ;
set INPUT_INDEL ;
for i in `seq 1 $CHUNK`; 
 do
  INPUT_SNP+="-I $OUTPUT/SNPs/filtered_snps.$CHR.$i.vcf.gz "; 
  INPUT_INDEL+="-I $OUTPUT/INDELs/filtered_indels.$CHR.$i.vcf.gz "; 
done

 ## For SNPs
   gatk GatherVcfs ${INPUT_SNP} -O $OUTPUT/SNPs_per_Chr/$CHR.SNPs.vcf.gz -R $REF ;
   tabix -f $OUTPUT/SNPs_per_Chr/$CHR.SNPs.vcf.gz ;
   gatk VariantsToTable -R $REF -V $OUTPUT/SNPs_per_Chr/$CHR.SNPs.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O $OUTPUT/SNPs_Var_Table/$CHR.table ;
 ## For INDELs
   gatk GatherVcfs ${INPUT_INDEL} -O $OUTPUT/INDELs_per_Chr/$CHR.INDELs.vcf.gz -R $REF ;
   tabix -f $OUTPUT/INDELs_per_Chr/$CHR.INDELs.vcf.gz ;
   gatk VariantsToTable -R $REF -V $OUTPUT/INDELs_per_Chr/$CHR.INDELs.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -O $OUTPUT/INDELs_Var_Table/$CHR.table ;
