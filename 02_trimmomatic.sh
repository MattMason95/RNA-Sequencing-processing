echo "Read Trimming {01_trimmomatic.sh}" # Declare start of trimming process

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
nTasks=$4   #Naming convention for raw reads

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
      -threads $nTasks \
      -phred33 \
      -basein $raw/${fq1} $raw/${fq2} \
      -baseout $trim/${output}_R1p.trimmed.fastq.gz $trim/${output}_R1u.trimmed.fastq.gz $trim/${output}_R2p.trimmed.fastq.gz $trim/${output}_R2u.trimmed.fastq.gz \
      ILLUMINACLIP:"$adapters":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

  # Option to remove the unpaired files
  # rm $trim/${output}_R1.U.fastq.gz \
  #    $trim/${output}_R2.U.fastq.gz

done

## LOGFILE
# sh ./AL_Logs.sh "log" "Trimmomatic" "$trim"

## FIN
