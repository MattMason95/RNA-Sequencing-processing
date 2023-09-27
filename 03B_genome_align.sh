#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Read Alignment with STAR {03_read_alignment.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

module star/2.5.0a

## PASSED VARIABLES 

trim=$1    #Location to read trimmed reads
bam=$2     #Location to write BAM files
fpkm=$3    #Location to write FPKM files
ctab=$4    #Location to write CTAB files
nTasks=$5  #Number of threads committed to the job

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ ALIGNMENT                            <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

for file1 in $trim/*_R1p.trimmed.fastq.gz 
do

file2=${file1%%_R1p.trimmed.fastq.gz}"_R2p.trimmed.fastq.gz" # Instantiate R2 singleton (PAIRED)

fq1=$(basename "${file1}")
fq2=$(basename "${file2}")
base=${fq1%%_R1p.trimmed.fastq.gz}

if find  "$fpkm" -name "*$base*"; then
  echo "$CHECKFILE exists. Skipping file."
else
        
  echo $fq1
  echo $fq2
  echo $base

  cmd2="STAR \
  --runThreadN $nTasks \
  --genomeDir $HOME/rds/hpc-work/Data/ENSEMBL/STAR \
  --readFilesIn ${fq1} ${fq2} \
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

  #echo $cmd2
  #eval $cmd2
fi 

done
wait

## FIN
