#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Pipeline module for STAR Genome indexing and read alignment

## SBATCH parameters
#SBATCH -p skylake-himem
#SBATCH -J mm2458_Read alignment
#SBATCH -A SPILLANTINI-SL3-CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=mm2458@cam.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=$HOME/rds/hpc-work/eofiles/%x.%j.out
#SBATCH --error=$HOME/rds/hpc-work/eofiles/%x.%j.err
#SBATCH --no-requeue

echo "Read Trimming {03_read_alignment.sh}" # Declare start of trimming process

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR

## Command line slurm submission script
bash Scripts/03_read_alignment.sh $HOME/rds/hpc-work/Data/Trimmed $HOME/rds/hpc-work/Data/BAM $HOME/rds/hpc-work/Data/FPKM $HOME/rds/hpc-work/Data/CTAB $SLURM_NTASKS
