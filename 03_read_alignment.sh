#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Read Alignment with STAR {03_read_alignment.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

module load $HOME/rds/hpc-work/Software/STAR
module load $HOME/rds/hpc-work/Software/samtools
module load $HOME/rds/hpc-work 

## PASSED VARIABLES

trim=$1    #Location to read trimmed reads
bam=$2     #Location to write BAM files
fpkm=$3    #Location to write FPKM files
ctab=$4    #Location to write CTAB files

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> GENOME BUILDING                           <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## FILE DOWNLOAD
cd $HOME/rds/hpc-work/Data
mkdir ENSEMBL
cd ENSEMBL
mkdir STAR 

## FTP FILE DOWNLOAD
## Use wget to download files for genome indexing
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #Primary Assembly of mouse genome
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz                       #Annotations of mouse genome
wget https://ftp.ensembl.org/pub/release-110/variation/vcf/mus_musculus/mus_musculus.vcf.gz                        #SNP variations of mouse genome

## UNZIP ENSEMBL FILES
gunzip -k Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip -k Mus_musculus.GRCm39.110.gtf.gz

## GENOME INDEXING
STAR --runThreadN 4 --runMoade genomeGenerate \
     --genomeDIR $HOME/rds/hpc-work/Data/ENSEMBL/STAR \
     --genomeSAindexNBases 10 \
     --sjdbGTFfile $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf \
     --genomeFastaFiles $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.dna.primary_assembly.fa

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ ALIGNMENT                            <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

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
