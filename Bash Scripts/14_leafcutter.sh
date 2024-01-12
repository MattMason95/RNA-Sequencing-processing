#!/bin/bash
# AUTHOR: Sarubini Kananathan, Matthew Mason
# EDITED: 01.2024
# DESC: Script for running leafcutter for alternative splicing analysis
echo "Performing alternative splicing analysis {14_leafcutter.sh}" # Declare start of trimming process

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> BUILD ENVIRONMENT                          <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

# Always add this command to your scripts
source $(conda info --base)/etc/profile.d/conda.sh

conda activate misopy

echo enviroment: $CONDA_DEFAULT_ENV

#$ -N leafcutter_DRG
#$ -o /lrlhps/users/c293525/Projects/Pain_ASO/src/logs/leafcutter_DRG.out
#$ -e /lrlhps/users/c293525/Projects/Pain_ASO/src/logs/leafcutter_DRG.err
#$ -cwd
#$ -l m_mem_free=220G
 
### Load packages
module load regtools/0.5.2
module load samtools/1.13
module load python/2.7.15
module load R/3.6.0

# <><><><><><><><><><><><><><><><><><><><><><><><>
# <> INSTANTIATE VARIABLES                      <>
# <><><><><><><><><><><><><><><><><><><><><><><><>

## Passed variables
SAMPLE="$1"

ALIGNED="$HOME/rds/hpc-work/Data/Aligned_reads"

for BAM in $ALIGNED/$SAMPLE; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $BAM
    regtools junctions extract -a 8 -s 2 -m 50 -M 500000 $BAM -o $BAM.junc
    echo $BAM.junc >> JunctionFiles.txt
done

### intron clustering (requires packages python/3.6.5 version)
python3 /lrlhps/users/c293525/Packages/leafcutter/clustering/leafcutter_cluster_regtools.py -j test_juncfiles.txt -m 50 -o p301sVSc57 -l 500000

 
 
## Differential intron excision analysis
# R --no-save -e 'install.packages("RcppParallel")'
# R --no-save -e 'install.packages("loo")'
# R --no-save -e 'install.packages("optparse")'
 
cd /lrlhps/data/neuroscience/common/leafcutter
./run_leafcutter.sh /opt/leafcutter/scripts/leafcutter_ds.R --num_threads 8 -e /lrlhps/data/neuroscience/common/leafcutter/Ensembl38_109_all_exons.txt.gz /lrlhps/genomics/prod/Neuroscience/Pain_Human_DRG_Lilly76/user_data/results/public_alignment/STAR/testPainvsControl_perind_numers.counts.gz /lrlhps/users/c293525/Projects/Pain_ASO/data/staging/groups_file.txt
mv leafcutter_ds_cluster_significance.txt leafcutter_ds_effect_sizes.txt /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter
# 
./run_leafcutter.sh /opt/leafcutter/scripts/ds_plots.R -e /lrlhps/data/neuroscience/common/leafcutter/Ensembl38_109_all_exons.txt.gz /lrlhps/genomics/prod/Neuroscience/Pain_Human_DRG_Lilly76/user_data/results/public_alignment/STAR/testPainvsControl_perind_numers.counts.gz /lrlhps/users/c293525/Projects/Pain_ASO/data/staging/groups_file.txt /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/leafcutter_ds_cluster_significance.txt -f 0.05
mv ds_plots.pdf /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter
 
# Changing gtf2toleafcutter (run only once)
# cd /lrlhps/data/neuroscience/common/leafcutter
#./run_leafcutter.sh /opt/leafcutter/leafviz/gtf2leafcutter.pl -o /lrlhps/data/neuroscience/common/leafcutter/Ensembl38_104 /lrlhps/users/c293525/Genomes/Homo_sapiens.GRCh38.104.mod.gtf
#./run_leafcutter.sh /opt/leafcutter/leafviz/gtf2leafcutter.pl -o /lrlhps/data/neuroscience/common/leafcutter/Ensembl38_109 /lrlhps/users/c293525/Genomes/Ensembl_v109/Homo_sapiens.GRCh38.109.gtf
 
 
./run_leafcutter.sh /lrlhps/users/c293525/Packages/leafcutter/leafviz/prepare_results.R -m /lrlhps/users/c293525/Projects/Pain_ASO/data/staging/groups_file.txt /lrlhps/genomics/prod/Neuroscience/Pain_Human_DRG_Lilly76/user_data/results/public_alignment/STAR/testPainvsControl_perind_numers.counts.gz /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/leafcutter_ds_cluster_significance.txt /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/leafcutter_ds_effect_sizes.txt Ensembl38_109 -o /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/ChronicPainvsControl.Rdata 
# ./run_leafcutter.sh /lrlhps/users/c293525/Packages/leafcutter/leafviz/prepare_results_ruby.R -m /lrlhps/users/c293525/Projects/Pain_ASO/data/staging/groups_file.txt /lrlhps/genomics/prod/Neuroscience/Pain_Human_DRG_Lilly76/user_data/results/public_alignment/STAR/testPainvsControl_perind_numers.counts.gz /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/leafcutter_ds_cluster_significance.txt /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/leafcutter_ds_effect_sizes.txt Ensembl38_109 -o /lrlhps/users/c293525/Projects/Pain_ASO/data/intermediate/leafcutter/ChronicPainvsControl.Rdata
 
 
#qsub -R y -pe smp 12 -M kananathan_sarubini@network.lilly.com -m be -cwd leafcutter.sh

