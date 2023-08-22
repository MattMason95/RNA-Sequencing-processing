#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Pipeline module for Trimmomatic trimming of reads

## SBATCH parameters
#SBATCH -p skylake-himem
#SBATCH -J mm2458_Trimmomatic
#SBATCH -A SPILLANTINI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH --mail-user=mm2458@cam.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=$HOME/rds/hpc-work/eofiles/%x.%j.out
#SBATCH --error=$HOME/rds/hpc-work/eofiles/%x.%j.err
#SBATCH --no-requeue

echo "Read Trimming {01_trimmomatic.sh}" # Declare start of trimming process

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> BUILD DEFAULT ENVIRONMENT                 <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

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

  output=${base}

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
  # rm $trim/${output}_R1.U.fastq.gz \
  #    $trim/${output}_R2.U.fastq.gz

done

## LOGFILE
# sh ./AL_Logs.sh "log" "Trimmomatic" "$trim"

## FIN
