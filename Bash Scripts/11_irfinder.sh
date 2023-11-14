#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Script for running irfinder for intron retention analysis

echo "Analyse intron rention {11_irfinder.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate bioinformatics

echo enviroment: $CONDA_DEFAULT_ENV

## PASSED VARIABLES 
sample="$1"

## STATIC VARIABLES 
REF="$HOME/rds/hpc-work/Data/ENSEMBL/


cd $HOME/rds/hpc-work/Data/Aligned_reads/$sample
mkdir -p irfinder

BAM=*sorted.bam.out

cmd1="irfinder -m BAM -r $REF -d irfinder $BAM"  

echo $cmd1
eval $cmd1
