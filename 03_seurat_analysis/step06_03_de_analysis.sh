#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/intermediate_data/log/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/intermediate_data/log/
#$ -pe smp 2
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

#make the results output directory
mkdir -p /gladstone/bioinformatics/adnetworksppg/Project_1/NK01/results/snRNA_mm10/analysis_hapoe_chr/05_seurat_analysis_filter_percent_mt_noGFAPCre/data/06_analysis_of_cluster_of_interest/enriched_genes/

data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/
script_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/scripts/snRNA_mm10/insert_HumanApoeTransgene_as_separate_chromosome
container_dir=/gladstone/bioinformatics/adnetworksppg/Project_1/NK01/parameters/containers
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/r_sn_rna_seq.sif Rscript $script_dir/seurat_analysis_noGFAPCre/src/06_03_de_analysis.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
