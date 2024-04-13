#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Genome Building with STAR {03A_genome_build.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate bioinformatics

echo enviroment: $CONDA_DEFAULT_ENV

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> VARIABLE DECLARATIONS                     <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

nTasks=24

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> FILE AND DIRECTORY PREPARATION            <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## FILE DOWNLOAD
cd Data
mkdir ENSEMBL
cd ENSEMBL
mkdir STAR 

## FTP FILE DOWNLOAD
## Use wget to download files for genome indexing
wget https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz #Primary Assembly of mouse genome
wget https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz                       #Annotations of mouse genome
wget https://ftp.ensembl.org/pub/release-110/variation/vcf/mus_musculus/mus_musculus.vcf.gz                        #SNP variations of mouse genome

## UNZIP ENSEMBL FILES
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.110.gtf.gz

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> GENOME BUILDING                           <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## GENOME INDEXING
cmd1="STAR --runThreadN $nTasks --runMode genomeGenerate \
--genomeDir $HOME/rds/hpc-work/Data/ENSEMBL/STAR \
--genomeFastaFiles $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.dna.primary_assembly.fa \
--sjdbGTFfile $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf \
--sjdbOverhang 149 \
--sjdbGTFtagExonParentTranscript Parent"

echo $cmd1
eval $cmd1

## Two-pass Genome Building
# cmd2="STAR --runThreadN $nTasks --runMode genomeGenerate \ 
# --genomeDir $HOME/rds/hpc-work/Data/ENSEMBL/Two_pass \ 
# --genomeFastaFiles $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.dna.primary_assembly.fa \
# --sjdbGTFfile $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf \
# --sjdbFileChrStartEnd SJ_out_filtered.tab \
# --sjdbOverhang 149"

# echo $cmd2
# eval $cmd2

## FIN
