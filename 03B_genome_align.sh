#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Read Alignment with STAR {03_read_alignment.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

module load star/2.5.0a

## PASSED VARIABLES 
fq1=$1
fq2=$2

## STATIC VARIABLES 
trim="$HOME/rds/hpc-work/Data/Trimmed"    #Location to read trimmed reads
bam="$HOME/rds/hpc-work/Data/BAM"     #Location to write BAM files
fpkm="$HOME/rds/hpc-work/Data/FPKM"    #Location to write FPKM files
ctab="$HOME/rds/hpc-work/Data/CTAB"    #Location to write CTAB files
nTasks=$5  #Number of threads committed to the job

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ ALIGNMENT                            <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

trimFQ1=${fq1%%_1.fq.gz}"_R1p.trimmed.fastq.gz"
trimFQ2=${fq2%%_2.fq.gz}"_R2p.trimmed.fastq.gz"
base=${fq1%%_1.fq.gz}

if find  "$fpkm" -name "*$base*"; then
  echo "$CHECKFILE exists. Skipping file."
else
        
  echo $trimFQ1
  echo $trimFQ2
  echo $base

  cmd2="STAR \
  --runThreadN 40 \
  --genomeDir $HOME/rds/hpc-work/Data/ENSEMBL/STAR \
  --readFilesIn $trim/$trimFQ1 $trim/$trimFQ2 \
  --readFilesCommand zcat \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix ${base}_ \
  --quantMode GeneCounts"

  echo $cmd2
  eval $cmd2
fi 

done
wait

## FIN
