#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -e /gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/tmp/log/
#$ -pe smp 2
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

data_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/
script_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/scripts/YH_NK01/insert_HumanApoeTransgene_as_separate_chromosome
container_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/assets/containers
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/r_sn_rna_seq.sif Rscript $script_dir/seurat_analysis_noGFAPCre/src/step14_01_de_gene_analyses_revision_reviewer2.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
