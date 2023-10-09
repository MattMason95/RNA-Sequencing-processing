echo "Read Trimming {02_trimmomatic.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ TRIMMING                             <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## PASSED VARIABLES
raw="$HOME/rds/hpc-work/Data/Fasta_files"                                        #Location to read raw reads
trim="$HOME/rds/hpc-work/Data/Trimmed"                                           #Location to write trimmed reads
adapters="$HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"   #Location to read adapaters   

## TRIMMOMATIC LOCATION
TRIMMOMATIC="$HOME/rds/hpc-work/Software/trimmomatic-0.39.jar"

## EXECUTION
fq1=$1
fq2=$2

base=${fq1%%_1.fq.gz}

CHECKFILE="${base}_R1u.trimmed.fastq.gz"

if [ -f $trim/$CHECKFILE ]; then
        echo "$CHECKFILE exists. Skipping file."
else
        echo $fq1
        echo $fq2
        echo $base

        ## Run trimmomatic
        java -jar "$TRIMMOMATIC" PE -threads 40 $raw/$fq1 $raw/$fq2 $trim/${base}_R1p.trimmed.fastq.gz $trim/${base}_R1u.trimmed.fastq.gz $trim/${base}_R2p.trimmed.fastq.gz $trim/${base}_R2u.trimmed.fastq.gz ILLUMINACLIP:"$adapters":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

wait
