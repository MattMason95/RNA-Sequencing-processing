#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Indexing of BAM files {04_bam_index.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate bioinformatics

echo enviroment: $CONDA_DEFAULT_ENV

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

echo sample: "$sample"

## FIND BAM FILE 
for bam in $ALIGNED/$sample/*.out.rd.rg.bam;
do

base=$(basename "$bam")
output=${$base%%.out.bam}".out.rd.rg.bam"

echo path: $bam
echo file: $base
echo root: $output

cd $ALIGNED/$sample

echo pwd: "$PWD"

## EXECUTE
cmd1="picard MarkDuplicates \
  I=$base \
  O=$output \
  M=$base.txt"

echo $cmd1
eval $cmd1
done
