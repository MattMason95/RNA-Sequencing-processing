#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Pipeline module for aligning trimmed reads to mouse genome index

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> SBATCH DIRECTIVES                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

#SBATCH -J mm2458_ReadAlignment
#SBATCH -A SPILLANTINI-SL3-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mm2458@cam.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --no-requeue

export nTasks=$SLURM_CPUS_PER_TASK

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> LOAD ENVIRONMENT                          <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

# ~~~~~~~~~~ DEFAULT REQUIREMENTS ~~~~~~~~~~
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-ccl              # REQUIRED - loads the basic environment

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> PREPARE EXECUTABLE                        <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

# ~~~~~~~~~~ EXECUTEABLE ~~~~~~~~~~
# Passed variables:
# 1 = Directory to READ trimmed reads
# 2 = Directory to WRITE BAM files 
# 3 = Directory to WRITE FPKM files
# 4 = Directory to WRITE CTAB files
# 5 = SLURM N-tasks
bash Scripts/03B_read_alignment.sh $HOME/rds/hpc-work/Data/Trimmed $HOME/rds/hpc-work/Data/BAM $HOME/rds/hpc-work/Data/FPKM $HOME/rds/hpc-work/Data/CTAB $nTasks

echo "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

## FIN
