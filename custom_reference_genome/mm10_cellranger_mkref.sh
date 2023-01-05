#!/bin/bash

#genome metadata
genome="adppg-mm10-apoe-chr-mapt-chr"
version="2021-A"

#paths
output_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters

#input files
fasta_custom="${output_dir}/adppg-mm10-apoe-chr-mapt-chr-2021-A_build/Mus_musculus.GRCm38.dna.primary_assembly.fa.custom"
gtf_custom="${output_dir}/adppg-mm10-apoe-chr-mapt-chr-2021-A_build/gencode.vM23.primary_assembly.annotation.gtf.custom"

cd $output_dir
cellranger mkref --ref-version="$version" --genome="$genome" --fasta="$fasta_custom" --genes="$gtf_custom"

