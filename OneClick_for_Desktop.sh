#!/bin/bash

## User environment 
export WORKDIR=`pwd` ;
export SIF=$WORKDIR/BioApps.sif ;
export CHRs=`cat $WORKDIR/ref/*.fasta.fai | wc -l` ;
COUNT=0;

## Check for singularity
if command -v singularity &> /dev/null 
then
 echo " Good! Singularity is available!" ;
 COUNT=$(( COUNT + 1 ));
else
 echo " Singularity is unavailable" ;
 echo " Your program is exiting!"  ;
 exit;
fi

## Check for MPI library 
if command -v mpicc &> /dev/null
then 
 echo " Excellent, MPI library available!"
 COUNT=$(( COUNT + 1 ));
else 
 echo " Please install MPI Library for parallel runs!"
 echo " Exiting your run now!"
fi
 
## Ensure both Singularity & MPI libraries are available 

if [ $COUNT == 2 ]; 
then
  export SINGULARITY_BIND="$WORKDIR,$PWD,/sw" ;
  ## Create required hostfile for mpi jobs
  NUM=`cat /proc/cpuinfo | grep processor | wc -l` ;
  for i in `seq 1 $NUM`;
   do
    echo `hostname` >> hostfile ;
   done
   echo " Host files are generated! and program will be executing now!" ;
  ## Phase 1
      cd $WORKDIR/input ;
      ls -lrta *.1.fastq.gz | awk '{print $9}' > $WORKDIR/scripts/Phase1/Phase1.txt ;
      SAMPLES=`cat $WORKDIR/scripts/Phase1/Phase1.txt | wc -l` ;
      cd $WORKDIR ;
      mpicc ${WORKDIR}/scripts/Phase1/Phase1.c -o ${WORKDIR}/scripts/Phase1/Phase1.exe
      mpiexec -np $SAMPLES -hostfile $WORKDIR/hostfile singularity exec $SIF ${WORKDIR}/scripts/Phase1/Phase1.exe

  ## Phase 2
      mpicc ${WORKDIR}/scripts/Phase2/Phase2.c -o ${WORKDIR}/scripts/Phase2/Phase2.exe ;
      mpiexec -np $SAMPLES -hostfile $WORKDIR/hostfile singularity exec $SIF ${WORKDIR}/scripts/Phase2/Phase2.exe ;

  ## Phase 3
     if [ $SAMPLES -lt $CHRs ]
     then 
        export SAMPLES=$CHRs ;
     fi 
     mpicc ${WORKDIR}/scripts/Phase3/Phase3.c -o ${WORKDIR}/scripts/Phase3/Phase3.exe ;
     sh ${WORKDIR}/scripts/Phase3/Phase3.prerequisite.optimized.sh $SAMPLES ;
     mpiexec -np $SAMPLES -hostfile $WORKDIR/hostfile singularity exec $SIF ${WORKDIR}/scripts/Phase3/Phase3.exe

  ## Phase 4
    mpicc ${WORKDIR}/scripts/Phase4/Phase4.c -o ${WORKDIR}/scripts/Phase4/Phase4.exe ;
    source ${WORKDIR}/scripts/Phase4/Phase4.prerequisite.sh ;
    mpiexec -np $SAMPLES -hostfile $WORKDIR/hostfile singularity exec $SIF ${WORKDIR}/scripts/Phase4/Phase4.exe ;

else
 echo " **** Please ensure Singularity and MPI libraries are available in your system ***" ;
fi
