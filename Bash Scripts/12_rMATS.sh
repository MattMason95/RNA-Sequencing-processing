#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 03.2024
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Performing Alternative Splicing Analysis with rMATS {15_rMATS.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate rMATS

echo $CONDA_DEFAULT_ENV

## STATIC VARIABLES 
GTF="/home/mm2458/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf"    #Location to read trimmed reads
TEMP="/home/mm2458/rds/hpc-work/Data/rMATS/rmats_tmp"      # Location for temporary files
OUTPUT="/home/mm2458/rds/hpc-work/Data/rMATS"     #Location to write output files

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## Allocate BAM txt files 
B1="/home/mm2458/rds/hpc-work/Data/rMATS/p301s_7m.txt"
B2="/home/mm2458/rds/hpc-work/Data/rMATS/c57_7m.txt"

## Execute rMATS
cmd1="~/miniconda3/envs/rMATS/bin/python ~/miniconda3/envs/rMATS/bin/rmats.py --gtf $GTF --b1 $B1 --b2 $B2 --od $OUTPUT --tmp $TEMP -t paired --readLength 150"

echo $cmd1
eval $cmd1

