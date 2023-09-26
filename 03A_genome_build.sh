#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Genome Building with STAR {03A_genome_build.sh}" # Declare start of trimming process

module star/2.5.0a

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> VARIABLE DECLARATIONS                     <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

nTasks=$1

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> FILE AND DIRECTORY PREPARATION            <>
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

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> GENOME BUILDING                           <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## GENOME INDEXING
cmd1 = "STAR --runThreadN $nTasks --runMode genomeGenerate \
--genomeDIR $HOME/rds/hpc-work/Data/ENSEMBL/STAR \
--genomeSAindexNBases 10 \
--sjdbGTFfile $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf \
--genomeFastaFiles $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.dna.primary_assembly.fa"

echo $cmd1
eval $cmd1 

## FIN
