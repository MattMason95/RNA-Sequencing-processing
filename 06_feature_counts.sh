#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 10.2023
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Counting of fragments per genomic locus {06_feature_counts.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate bioinformatics

echo enviroment: $CONDA_DEFAULT_ENV

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

echo sample: "$sample"

## FIND BAM FILE 
for bam in $ALIGNED/$sample/*.out.rd.rg.bam;
do

base=$(basename "$bam")
output=${$base%%.out.bam}".out.rd.rg.bam"

echo path: $bam
echo file: $base
echo root: $output

cd $ALIGNED/$sample

echo pwd: "$PWD"

## EXECUTE FEATURECOUNTS
cmd1="featureCounts \
  -p --countReadPairs \
  -a $HOME/rds/hpc-work/Data/ENSEMBL/Mus_musculus.GRCm39.110.gtf \
  -t exon \
  -g gene_id \
  -o ${base}_featureCounts.txt \
  $base"

echo $cmd1
eval $cmd1

done
