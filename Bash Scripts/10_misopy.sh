#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Script for running misopy in parallel on the cluster for alternative splicing analysis

echo "Performing alternative splicing analysis {10_misopy.sh}" # Declare start of trimming process

## Settings file location: /home/mm2458/miniconda3/envs/misopy/lib/python2.7/site-packages/misopy/settings

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate misopy

echo enviroment: $CONDA_DEFAULT_ENV

## PASSED VARIABLES 
sample="$1"

## STATIC VARIABLES 
ALIGNED="$HOME/rds/hpc-work/Data/Aligned_reads"    #Location to read trimmed reads
BAI="$HOME/rds/hpc-work/Data/BAI"     #Location to write BAM files

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> CONVERT GTF to GFF3                        <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## FILETYPE CONVERSION
GTF="$HOME/rds/hpc-work/Data/ENSEMBL/"
GFF3="Mus_musculus.GRCm39.110.gff3"

if [ -f $GTF/$GFF3 ]; then
   echo "GFF3 file - $GFF3 - already exists. Proceeding with MISO"
else
   echo "GFF3 file does not exist. Converting GTF to GFF3"
   
   ## Utilise AGAT for conversion
   
   cd $GTF
   cmd1="agat_convert_sp_gxf2gxf.pl -g Mus_musculus.GRCm39.110.gtf -o $GFF3"
   echo $cmd1
   #eval $cmd1
   wait
   [ -f $GFF3 ] && echo "GFF3 file created!" 
fi

## GFF INDEXING
mkdir -p indexed_gff

cmd2="miso index_gff --index $GFF3 indexed_gff"
echo $cmd2 
#eval $cmd2
wait

cd ../..
echo "$PWD"

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> GFF INDEXING                               <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

cd Data
mkdir -p miso_output
cd miso_output

BAM="$HOME/rds/hpc-work/Data/Aligned_reads/$sample/*.out.rd.bam"
IndexDIR="$HOME/rds/hpc-work/Data/ENSEMBL/indexed_gff"

cmd3="miso --run $IndexDIR $BAM --output-dir ./miso_output --read-len 150 --paired-end 250 15 --use-cluster"

echo $cmd3 
#eval $cmd3
done

## FIN 
