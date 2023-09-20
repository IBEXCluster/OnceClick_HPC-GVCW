#!/bin/bash
export REF=${WORKDIR}/ref/Nipponbare_chr.fasta ;
export INPDIR=${WORKDIR}/output/VCF ;
export MaxRanks=$1 ;
#export WORKDIR=$PWD
function main()
{
#  read -p " How many CPUs (cores) you are planning to use:  " MaxRanks
#  read -p " Please enter absolute path of the reference genome with name (e.g. /project/k01/kathirn/3k/ref/Nipponbare_chr.fasta) : " REF
#  read -p " Please enter absolute path of your gVCF files directory: " INPDIR 
  Validate;
  FindMaxSizeofChr ;
  FindOptimalChunk ;
  DistributeParallel ;
  TasksPerNode ;
  InputPrepare ;
}

function Validate()
{
  ## Check the number of CPUs are integer 
  if ! [[ "$MaxRanks" =~ ^[0-9]+$ ]]; then
   echo "Please give number of CPUs (in number ONLY)"
   exit;    
  fi

  ## Check the Reference Genome Index
  if [ ! -f $REF.fai ]; then
    echo "Reference Genome Index file not found! Please create Index file ;-)"
    exit 1;
  fi
}

function FindMaxSizeofChr()
{
  MaxSize=0 ;
  while IFS=$'\t' read -r -a myREF
  do
     ChrName=${myREF[0]};
     ChrLen=${myREF[1]};
     if [ $MaxSize -lt $ChrLen ]; then
       MaxSize=$ChrLen
     fi
  done < $REF.fai 
   echo " "
   echo "*************************************************************"
  echo "Max size of Chromosome in the given reference is:" $MaxSize
  #Find the total numbers of Chr
  TotalChr=`cat $REF.fai | wc -l `
  echo "Total no. of Chromosomes in the given reference is:"  $TotalChr 
  if [[ $MaxRanks -lt $TotalChr ]]; then
   echo "This algorithm required at least $TotalChr CPUs. However, the provided CPUs=$MaxRanks is insufficient for performance improvement"
    exit 1;
   echo "*************************************************************"
  else
   echo "No. of optimal CPUs will be calculated as follows:"
#   echo " ... Using 4 Cores/Sample and 8 samples/Node ... (in Intel Sakylake architecture) ..."
  fi
}


function FindOptimalChunk()
{
  #Find the best possible Chunk size 
  Divide=1;
  AssinedJobs=1;
  while [ ${MaxRanks} -gt $AssinedJobs ];
  do
     Divide=$(( $Divide + 1 ));
     AssinedJobs=$(( $Divide * $TotalChr )) ;
  done
  Divide=$(( Divide - 1 ));
  Chunk=$(( $MaxSize / $Divide ));
  #echo "total Chr = $TotalChr, Max Chr length = $MaxSize, Chunk number = $Divide, Chunk size = $Chunk "
}



function DistributeParallel()
{ 
  #Prepare Chromosome Chunk intervals in a file "submitted_chunks_data.txt"
  rm -Rf ${WORKDIR}/scripts/Phase3/submitted_chunks_data.txt ${WORKDIR}/scripts/Phase3/split.txt ${WORKDIR}/scripts/Phase3/distribution.txt ${WORKDIR}/scripts/Phase3/Phase3.input.list;
  while IFS=$'\t' read -r -a myREF
  do
    ChrName=${myREF[0]};
    ChrLen=${myREF[1]};
    Part=$(( $ChrLen / $Chunk )) ;
    tmp=$(( $Part * $Chunk )) ;
    if [ $ChrLen -gt $tmp ]; then
      Part=$(( $Part + 1 ));
      ## Equal partitiaons 
      echo "$ChrName split into $Part parts" >> ${WORKDIR}/scripts/Phase3/submitted_chunks_data.txt
    else
      ## Last part of the partition
      echo "$ChrName split into $Part parts" >> ${WORKDIR}/scripts/Phase3/submitted_chunks_data.txt
    fi
   done < $REF.fai

  #Prepare a Chromosome Interval Distributoion across all the reference genome 
   while IFS=$'\t' read -r -a myREF
   do
     size=1;
     DefSize=$Chunk;
     ChunkSize=$DefSize;
     ChrName=${myREF[0]};
     ChrLen=${myREF[1]};
  
     ## PREPARE CHUNKS [Begin:End] Intervals for EACH CHROMOSOME 
      for ((Start=1; Start<$ChrLen; Start+=$DefSize))
      do
         End=$(( $size*$ChunkSize ))
         if [ $End -lt $ChrLen ]
         then
           echo "$ChrName   $size   $Start  $End"  >> ${WORKDIR}/scripts/Phase3/split.txt
           size=$(( $size + 1 ));
         else   ### Last chunk 
           Start=$(( $End - $ChunkSize +1 ));
           End=$ChrLen;
           echo "$ChrName  $size  $Start  $End"   >> ${WORKDIR}/scripts/Phase3/split.txt
         fi ## END OF ALL CHUNKS within the specific CHROMOSOME 
          #echo "Preparing for the next Chunk:" $ChrName.chunk_$size 
       done ## END OF CHROMOSOME (Chr by Chr) 
       echo "$ChrName    $size" >> ${WORKDIR}/scripts/Phase3/distribution.txt
    done < $REF.fai
    ## END OF ALL CHROMOSOMEs
}

function TasksPerNode()
{ 
  TotalLines=`cat ${WORKDIR}/scripts/Phase3/split.txt | wc -l` ;
  #Nodes=$(( ((TotalLines+31) / 32) ))
#  Nodes=$(( ((TotalLines+31) / 4) ))
#  echo "The following values can be set in your Phase3.batch "
#  echo " -N $Nodes and --ntasks=$TotalLines " 
   echo "Please use: -np $TotalLines"
   echo "*************************************************************"
}

function InputPrepare()
{
 ## Add all the *.g.vcf files for CombineGVCF
 set INPUT_TMP
 for i in `ls -l ${INPDIR}/*.g.vcf.gz | awk '{print $9}'`
 do
   INPUT_TMP+="$i -V ";
 done
 INPUT=${INPUT_TMP::-4}
  echo $INPUT > ${WORKDIR}/scripts/Phase3/Phase3.input.list
}

main "$@"
