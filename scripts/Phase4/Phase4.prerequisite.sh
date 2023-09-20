#!/bin/bash
if [ ! -f ${WORKDIR}/scripts/Phase3/distribution.txt ]; then
    echo "distribution.txt ...File not found!"
    echo "We can't combine the Chunks of SNPs and INDELs" 
    break;
else
   cp ${WORKDIR}/scripts/Phase3/distribution.txt ${WORKDIR}/scripts/Phase4/distribution.txt ;
   count=`cat ${WORKDIR}/scripts/Phase4/distribution.txt | wc -l` ;
   echo "The number of processes for data (or parallel) distriibution is: -np $count" ;
   echo " Please use: -np $count in the mpirun" ;
fi 
