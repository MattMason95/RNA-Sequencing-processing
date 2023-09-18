#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Pipeline module for Trimmomatic trimming of reads

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> SBATCH DIRECTIVES                         <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

#SBATCH -J mm2458_Trimmomatic
#SBATCH -A SPILLANTINI-SL3-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -D $HOME/rds/hpc-work 
#SBATCH --time=02:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=mm2458@cam.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=eofiles/02_trimmomatic.log
#SBATCH --no-requeue

# ~~~~~~~~~~ COMPUTE ALLOCATIONS  ~~~~~~~~~~
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> LOAD ENVIRONMENT                          <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

# ~~~~~~~~~~ DEFAULT REQUIREMENTS ~~~~~~~~~~
#. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
#module purge                               # Removes all modules still loaded
#module load rhel8/default-icl              # REQUIRED - loads the basic environment

#module load pkg-config-0.29.2-gcc-6.2.0-we4glmw
#module load gcc/9

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> PREPARE EXECUTABLE                        <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

# ~~~~~~~~~~ EXECUTEABLE ~~~~~~~~~~
APP="Scripts/02_trimmomatic.sh"

ARGS="$HOME/rds/hpc-work/Data/Fasta_files $HOME/rds/hpc-work/Data/Trimmed $HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa $SLURM_NTASKS"

workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

CMD="$APP $ARGS"

export OMP_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> EXECUTE SBATCH                            <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
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

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD

## FIN
