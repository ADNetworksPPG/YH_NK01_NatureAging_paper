#!/bin/bash
#$ -cwd
#$ -o ~/log/
#$ -e ~/log/
#$ -pe smp 4
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=08:00:00
#$ -j yes

data_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr
container_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/containers
export SINGULARITY_BINDPATH="$data_dir"

cd $data_dir/results/01_cellranger_count/
singularity exec $container_dir/cellranger.sif cellranger count \
--id=S1_PS19-fE4_GFAP-Cre32 \
--transcriptome=$data_dir/parameters/adppg-mm10-apoe-chr-mapt-chr \
--fastqs=$data_dir/data/fastq/ \
--sample=PS19-fE4_GFAP-Cre32 \
--include-introns=true \

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
