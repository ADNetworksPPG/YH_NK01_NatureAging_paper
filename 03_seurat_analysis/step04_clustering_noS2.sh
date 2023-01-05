#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/intermediate_data/log/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/intermediate_data/log/
#$ -pe smp 1
#$ -l mem_free=80G
#$ -l scratch=50G
#$ -l h_rt=08:00:00
#$ -j yes

#make the results output directory
mkdir -p /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/results/snRNA_mm10/analysis_hapoe_chr/05_seurat_analysis_filter_percent_mt_noGFAPCre/data/04_clustering_noS2/
mkdir -p /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/results/snRNA_mm10/analysis_hapoe_chr/05_seurat_analysis_filter_percent_mt_noGFAPCre/plot/04_clustering_noS2/

data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/
script_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/scripts/snRNA_mm10/insert_HumanApoeTransgene_as_separate_chromosome
container_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/parameters/containers
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/r_sn_rna_seq.sif Rscript $script_dir/seurat_analysis_noGFAPCre/src/04_clustering_noS2.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
