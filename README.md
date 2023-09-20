# One Click HPC-GVCW 

## HPC-based genome variant calling workflow (HPC-GVCW)

Yong Zhou, Nagarajan Kathiresan, Zhichao Yu, Luis F. Rivera, Manjula Thimma, Keerthana Manickam, Dmytro Chebotarov, Ramil Mauleon, Kapeel Chougule, Sharon Wei, Tingting Gao, Carl D. Green, Andrea Zuccolo, Doreen Ware, Jianwei Zhang, Kenneth L. McNally, Rod A. Wing

bioRxiv 2023.06.25.546420; </br>

doi: https://doi.org/10.1101/2023.06.25.546420

## Overview 

 Here, we are providing One Click workflow for (a) any High-Performance Computing (HPC) plarform and (b) High-end workstations (HW).
 Using this one click, the 4 Phases of the rice genome workflow (https://github.com/IBEXCluster/Rice-Variant-Calling) will be executed automatically without manual processing. 
 
## Prerequisites 
 Singularity for executing our bio-container software stack
 MPI Library for parallel data processing across multiple-cores across the HPC or HW platform 
 Slurm scheduler for HPC platform. 

## Steps: 

**1. Clone the repository into your HPC or HW platform**
   
    git clone --recursive https://github.com/IBEXCluster/OnceClick_HPC-GVCW.git
    

**2. Change to ROOT directory**

    cd OneClick_HPC-GVCW/
    
**3. Singularity image file will be created using the following command:**
   
    singularity build BioApps.sif docker://ibexcluster/bioapps:v1.0
   
**4. Prepare Input and Reference files**

     The input and reference files will be prepared in the corresponding directories input/ and ref/ respectively. 
     More informations are availble at input/README and ref/README for better understanding. 
     

**5. Execute any one of the script:**

     For HPC platform: sbatch OneClick_for_HPC.sh
     For high-end workstation: sh ./OneClick_for_Desktop.sh
     

# Disclaimer 

This version is our preliminary **OneClick workflow** for testing. 
In case of any issues, please contact me at nagarajan.kathiresan@kaust.edu.sa 

 
