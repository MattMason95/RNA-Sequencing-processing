#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Indexing of BAM files {04_bam_index.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

module load samtools/1.10 

## PASSED VARIABLES 
sample=$1

## STATIC VARIABLES 
ALIGNED="$HOME/rds/hpc-work/Data/Aligned_reads"    #Location to read trimmed reads
BAI="$HOME/rds/hpc-work/Data/BAI"     #Location to write BAM files

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## ACCESS ARRAY-SPECIFIC SAMPLE DIR. 
cd $ALIGNED/$sample

## FIND BAM FILE 
bam=$*.out.bam

echo $bam

## EXECUTE
cmd1="samtools index $bam --output $BAI"

echo $cmd1
eval $cmd1

## FIN
