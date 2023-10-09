#!/bin.bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Read quality control with FASTQC {01_fastqc.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

module load fastqc/0.11.4

## PASSED VARIABLES

trim="$HOME/rds/hpc-work/Data/Trimmed/Paired"    #Location to read trimmed reads

cd $trim
#mkdir FastQC

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ QUALITY CONTROL                       <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

cmd1="fastqc -t 32 *.trimmed.fastq.gz"

echo $cmd1
eval $cmd1

mv *_fastqc* FastQC

## FIN
