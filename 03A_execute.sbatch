#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Pipeline module for building mouse genome index with STAR
# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> SBATCH DIRECTIVES                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

#SBATCH -J mm2458_GenomeBuild
#SBATCH -A SPILLANTINI-SL3-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=06:00:00
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
bash Scripts/03A_genome_build.sh $nTasks

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
