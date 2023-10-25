## RNA-Sequencing-processing
Consolidation of BASH and R scripts for the pre-processing of raw RNA-sequencing data.

The Bash Script directory contains an array of scripts for performing the pre-processing and alignment of raw FastQ files. 
Each script has an accompanying "*_execute.sbatch" script for deploying the scripts to the Cambridge CSD3 using SLURM Job Scheduling.

## Libraries used: 
FastQC -> Read Quality Control
STAR -> Genome Indexing; Read Alignment
SamTools -> Coverage Statistics
Picard -> Duplicate Removal
FeatureCounts -> Gene Read Counts

The R Script directory contains the scripts used for analysing the read count data generated by FeatureCounts, and combined using the "combine_featurecounts.R".

## Libraries used:
DeSeq2 -> Read Normalisation and Diff. Exp. Analysis
