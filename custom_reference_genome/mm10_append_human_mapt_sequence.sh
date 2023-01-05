#!/bin/bash

# Set up source and build directories
build="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/adppg-mm10-apoe-chr-mapt-chr-2021-A_build"

# set paths for input files
fasta_transgene="${build}/Mus_musculus.GRCm38.dna.primary_assembly.fa.transgene"
fasta_custom="${build}/Mus_musculus.GRCm38.dna.primary_assembly.fa.custom"
fasta_mapt="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/Human_MAPT/Human_MAPT.fasta"
gtf_transgene="${build}/gencode.vM23.primary_assembly.annotation.gtf.transgene"
gtf_custom="${build}/gencode.vM23.primary_assembly.annotation.gtf.custom"
gtf_mapt="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/Human_MAPT/Human_MAPT.gtf"

#append human MAPT fasta sequence at the end of the reference fasta file
cp $fasta_transgene $fasta_custom 
cat $fasta_mapt >> $fasta_custom
grep ">" $fasta_custom

#append human MAPT GTF file at the end of the reference GTF file
cp $gtf_transgene $gtf_custom
cat $gtf_mapt >> $gtf_custom
tail $gtf_custom


################### END ###################
