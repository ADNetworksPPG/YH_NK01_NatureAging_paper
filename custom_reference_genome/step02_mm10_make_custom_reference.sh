#!/bin/bash
#$ -cwd
#$ -pe smp 2
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=02:00:00
#$ -j yes

data_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/
script_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/scripts/custom_reference_genome
container_dir=/wynton/group/gladstone/biocore/projects/ADPPG/Huang/NicoleKoutsodendris-NK01-snRNAseq-pathology-mm10-oct-2021/analysis_hapoe_chr/parameters/containers
export SINGULARITY_BINDPATH="$data_dir"


#modify and filter source file 
$script_dir/mm10_modify_source_files.sh
echo "*********   modify done!  *************"

#add the custom human apoE to the modified source files
if [ "$?" -eq 0 ]; then
    #append the human apoe transgene sequence
    echo "*********   human apoE append starting!  *************"
    $script_dir/mm10_append_human_apoe_transgene_sequence.sh
    echo "*********   human apoE append done!  *************"
    if [ "$?" -eq 0 ]; then
    	#append the human MAPT sequence
        echo "*********  human MAPT append starting!  *************"
    	$script_dir/mm10_append_human_mapt_sequence.sh
        echo "*********  human MAPT append done!  *************"
    	if [ "$?" -eq 0 ]; then
    		#run the cellranger mkref to generate the custom reference genome
                singularity exec $container_dir/cellranger.sif $script_dir/mm10_cellranger_mkref.sh
                echo "************* end  *************!"
        fi
    fi
fi


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
