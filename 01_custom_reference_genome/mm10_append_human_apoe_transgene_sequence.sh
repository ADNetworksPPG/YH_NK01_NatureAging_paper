#!/bin/bash

# Set up source and build directories
build="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/adppg-mm10-apoe-chr-mapt-chr-2021-A_build"

# set paths for input files
fasta_modified="${build}/Mus_musculus.GRCm38.dna.primary_assembly.fa.modified"
fasta_transgene="${build}/Mus_musculus.GRCm38.dna.primary_assembly.fa.transgene"
fasta_hapoe="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/Human_apoE/hapoE_transgene.fa"
gtf_filtered="${build}/gencode.vM23.primary_assembly.annotation.gtf.filtered"
gtf_transgene="${build}/gencode.vM23.primary_assembly.annotation.gtf.transgene"
gtf_hapoe="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/Human_apoE/hapoE_transgene.gtf"

#append human MAPT fasta sequence at the end of the reference fasta file
cp $fasta_modified $fasta_transgene 
cat $fasta_hapoe >> $fasta_transgene
grep ">" $fasta_transgene

#append human MAPT GTF file at the end of the reference GTF file
cp $gtf_filtered $gtf_transgene
cat $gtf_hapoe >> $gtf_transgene
tail $gtf_transgene


################### END ###################
