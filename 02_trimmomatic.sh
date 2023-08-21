#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Read Trimming {01_trimmomatic.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ TRIMMING                             <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## PASSED VARIABLES

raw=$1         #Location to read raw reads
trim=$2        #Location to write trimmed reads
adapters=$3    #Location to read adapaters
# trimNAMES=$4   #Naming convention for raw reads 

## TRIMMOMATIC LOCATION

TRIMMOMATIC="$HOME/rds/hpc-work/Software/trimmomatic-0.39.jar"

## EXECUTION

for infile in $raw/*_1.fq.gz 
do

  base=$(basename ${infile} _1.fq.gz) # Get file basename
  fq1=${base}_1.fq.gz # Instantiate R1 singleton
  fq2=${base}_2.fq.gz # Instantiate R2 singleton (PAIRED)

  output=${base_}

  # Run trimmomatic
  java -jar "$TRIMMOMATIC" \
      PE \
      -threads 16 \
      -phred33 \
      ${fq1} \
      ${fq2} \
      $trim/${output}_R1.P.fastq.gz \
      $trim/${output}_R1.U.fastq.gz \
      $trim/${output}_R2.P.fastq.gz \
      $trim/${output}_R2.P.fastq.gz \
      ILLUMINACLIP:"$adapters":2:30:10 \
      LEADING:3 \
      TRAILING:3 \
      SLIDINGWINDOW:4:15 \
      MINLEN:36

  # Option to remove the unpaired files
  rm $trim/${output}_R1.U.fastq.gz \
     $trim/${output}_R2.U.fastq.gz


done

## LOGFILE
sh ./AL_Logs.sh "log" "Trimmomatic" "$trim"

## FIN
