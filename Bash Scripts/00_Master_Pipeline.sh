#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of modular BASH scripts for RNA-seq read alignment and quantification

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD MASTER ENVIRONMENT                  <>
# <><><>><><><><><><><><><><><><><><><><><><><><>
## Installation and module import requirements: 
## - Trimmomatic 
## - STAR
## - 
## - 

# ~~~~~~~~~~ DEFAULT REQUIREMENTS ~~~~~~~~~~
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

# ~~~~~~~~~~ PIPELINE SPECIFIC ~~~~~~~~~~
module load $HOME/rds/hpc-work/Software/STAR
module load $HOME/rds/hpc-work/Software

TRIMMOMATIC="$HOME/rds/hpc-work/Software/trimmomatic-0.39.jar"

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> INPUT VARIABLE DECLARATION                <>
# <><><>><><><><><><><><><><><><><><><><><><><><>
## A space to instantiate additional input arguments from the .SBATCH EXECUTE script 

nTasks=$1 # Number of threads allocated to the entire pipeline at submission

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> DIRECTORY ALLOCATION                      <>
# <><><>><><><><><><><><><><><><><><><><><><><><>
## Modifiable directory allocation for data read/write

# ~~~~~~~~~~ GENERAL ~~~~~~~~~~
raw="$HOME/rds/hpc-work/Data/Fasta_files"  # Location of the raw Fastq files

# ~~~~~~~~~~ TRIMMOMATIC ~~~~~~~~~~
adapters="$HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa" # Location of Illumina adapter sequences for clipping
trimmed="$HOME/rds/hpc-work/Data/Trimmed"                                      # Location to write/read trimmed Fastq files 

# ~~~~~~~~~~ STAR ~~~~~~~~~~
genome="$HOME/rds/hpc-work/Data/Genome"            # Location of the importred genome build for indexing
                                                   # At the time of running, Mus musculus genome build: ____
genomeIndex="$HOME/rds/hpc-work/Data/Genome/Index" # Location of the indexed genome
bam="$HOME/rds/hpc-work/Data/BAM"                  #Location to write BAM files
fpkm="$HOME/rds/hpc-work/Data/FPKM"                #Location to write FPKM files
ctab="$HOME/rds/hpc-work/Data/CTAB"                #Location to write CTAB files

# ~~~~~~~~~~ FASTQC ~~~~~~~~~~
QC="$HOME/rds/hpc-work/Data/FastQC"

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BEGIN MODULAR PIPLINE                     <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

# ~~~~~~~~~~ FASTQC ~~~~~~~~~~
echo "Initiating FastQC"
cmdQC="$HOME/rds/hpc-work/Software/FastQC/..."
echo $cmdQC
eval $cmdQC

# ~~~~~~~~~~ TRIMMOMATIC ~~~~~~~~~~
echo "Initiating Trimmomatic"
adapters="$HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"

cmd1="Scripts/02_trimmomatic.sh $raw $trimmed $adapters $nTasks"
echo $cmd1 
eval $cmd1 

# ~~~~~~~~~~ FASTQC ~~~~~~~~~~
echo "Initiating FastQC"
cmdQC="$HOME/rds/hpc-work/Software/FastQC/..."
echo $cmdQC
eval $cmdQC

# ~~~~~~~~~~ STAR ~~~~~~~~~~
echo "Initiating STAR Alignment"
cmd2="Scripts/03_read_alignment.sh $trimmed $trimmed $adapters $nTasks"
echo $cmd2 
eval $cmd2 

# ~~~~~~~~~~ STAR ~~~~~~~~~~
echo "Initiating STAR Alignment"
cmd2="Scripts/03_read_alignment.sh $trimmed $trimmed $adapters $nTasks"
echo $cmd2 
eval $cmd2 



echo "Initiating ${}"
echo "Initiating ${}"
echo "Initiating ${}"
echo "Initiating ${}"
echo "Initiating ${}"
