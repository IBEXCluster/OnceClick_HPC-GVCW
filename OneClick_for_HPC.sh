#!/bin/bash
#SBATCH --ntasks=10
#SBATCH -N 1
#SBATCH --mem-per-cpu=16GB
#SBATCH -J Singularity
#SBATCH --error=Singularity.log.%J.err
#SBATCH --output=Singularity.log.%J.out
#SBATCH --time=2-1:00:00

## User environment 
export WORKDIR=`pwd` ;
export SIF=$WORKDIR/BioApps.sif ;
export CHRs=`cat $WORKDIR/ref/*.fasta.fai | wc -l` ;

## Module file 
module load mpich/3.3/gnu-6.4.0 singularity 
export SINGULARITY_BIND="$WORKDIR,$PWD,/sw" ;

## Create required files and directories  
rm -Rf $WORKDIR/hostfile ;
scontrol show hostnames > $WORKDIR/hostfile ;

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
