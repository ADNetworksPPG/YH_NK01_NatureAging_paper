#this script should be run on the data transfer (dt) node on wynton
#!/bin/bash

# Create the reference_sources folder for the source files
source="/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/reference_sources"
mkdir -p "$source"

#specify the source files to be downloaded
fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
fasta_download="${source}/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz"
gtf_download="${source}/gencode.vM23.primary_assembly.annotation.gtf.gz"

# Download source files if they do not exist in reference_sources folder
cd $source
if [ ! -f "$fasta_download" ]; then
    curl -sSO "$fasta_url"
    if [ "$?" -eq 0 ]; then
        echo "Fasta file download completed!"
    else
        echo "ERROR: Fasta file download failed."
    fi
else
	echo "Fasta file already exists."
fi
if [ ! -f "$gtf_download" ]; then
    curl -sSO "$gtf_url"
    if [ "$?" -eq 0 ]; then
        echo "GTF file download completed!"
    else
        echo "ERROR: GTF file download failed."
    fi
else
	echo "GTF file already exists."
fi


################### END ###################
