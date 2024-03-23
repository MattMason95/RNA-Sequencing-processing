#!/bin/bash
# AUTHOR: Matthew Mason
# EDITED: 03.2024
# DESC: Full assembly pipeline of BASH scripts for RNA-seq read alignment and quantification

echo "Generating splicing Sashimi plots with rMATS2SASHIMI {16_sashimi.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate rMATS

echo $CONDA_DEFAULT_ENV

## PASSED VARIABLES
EVENT=$1
FILE=$2

echo $EVENT
echo $FILE

## STATIC VARIABLES
eFILE="/home/mm2458/rds/hpc-work/Data/rMATS/7m_sashimi/Sig_rMATS/$FILE"    #Location to read trimmed reads
TEMP="/home/mm2458/rds/hpc-work/Data/rMATS/7m__sashimi"      # Location for temporary files
OUTPUT="/home/mm2458/rds/hpc-work/Data/rMATS/7m_sashimi/$EVENT"     #Location to write output files

echo $eFILE
echo $OUTPUT

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> READ INDEXING                              <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## Allocate BAM txt files
B1="/home/mm2458/rds/hpc-work/Data/rMATS/7m_sashimi/p301s_7m.txt"
B2="/home/mm2458/rds/hpc-work/Data/rMATS/7m_sashimi/c57_7m.txt"
GROUP='/home/mm2458/rds/hpc-work/Data/rMATS/7m_sashimi/7m_grouping.gf'

## Execute rMATS
cmd1="rmats2sashimiplot --b1 $B1 --b2 $B2 --event-type $EVENT -e $eFILE --l1 7m_P301S --l2 7m_C57 --exon_s 1 --intron_s 5 -o $OUTPUT --color '#A88FAC,#468189' --group-info $GROUP"

echo $cmd1
eval $cmd1


