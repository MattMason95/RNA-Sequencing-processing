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

## Command line slurm submission script
bash Scripts/02_trimmomatic.sh $HOME/rds/hpc-work/Data/Fasta_files $HOME/rds/hpc-work/Data/Trimmed $HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa $SLURM_NTASKS
