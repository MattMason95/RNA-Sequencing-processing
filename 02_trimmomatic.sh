echo "Read Trimming {02_trimmomatic.sh}" # Declare start of trimming process

# <><><>><><><><><><><><><><><><><><><><><><><><>
# <> READ TRIMMING                             <>
# <><><>><><><><><><><><><><><><><><><><><><><><>

## PASSED VARIABLES
raw="$HOME/rds/hpc-work/Data/Fasta_files"                                        #Location to read raw reads
trim="$HOME/rds/hpc-work/Data/Trimmed"                                           #Location to write trimmed reads
adapters="$HOME/rds/hpc-work/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"   #Location to read adapaters
nTasks=$1                                                                        #Naming convention for raw reads

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
done
wait

## FIN
