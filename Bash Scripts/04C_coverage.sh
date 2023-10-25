#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Generate Coverage Stats of BAM files {04C_Coverage.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate bioinformatics

echo $CONDA_DEFAULT_ENV

## PASSED VARIABLES 
sample="$1"

## STATIC VARIABLES 
ALIGNED="$HOME/rds/hpc-work/Data/Aligned_reads"    #Location to read trimmed reads
BAI="$HOME/rds/hpc-work/Data/BAI"     #Location to write BAM files

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

echo sample: "$sample"

## FIND BAM FILE 
for bam in $ALIGNED/$sample/*.out.bam;
do

base=$(basename "$bam")

echo $bam
echo $base

cd $ALIGNED/$sample

echo "$PWD"

## EXECUTE
cmd1="samtools coverage $base"

echo $cmd1
eval $cmd1
done
