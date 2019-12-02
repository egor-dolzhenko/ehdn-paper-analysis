#!/bin/bash

# PBS Directives, assuming PBSpro
#PBS -A st-wasserww-1
#PBS -V
#PBS -N STR_EHDN_MiniSimFiles
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=24:0:0

EHDN8=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo-v0.8.0-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.0
EHDN62=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo-v0.6.2
GENOME_FASTA=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.fa
STR_ANALYSIS_DIR=/scratch/st-wasserww-1/STR_Analysis/
STR_BAM_DIR=$STR_ANALYSIS_DIR/STR_Simulation/str_bams/GRCh37_20191030/
STR_JSON_DIR=$STR_ANALYSIS_DIR/STR_Simulation/str_jsons/GRCh37_20191030/
mkdir -p $STR_JSON_DIR
cd $STR_JSON_DIR

for SAMPLE_BAM in $(ls $STR_BAM_DIR/*.sorted.bam )
do	

        echo $SAMPLE_BAM
        IFS='/' read -a array <<< $SAMPLE_BAM
        SampleBam=${array[-1]}
        IFS='.' read -a array2 <<< "$SampleBam"
        SAMPLE_ID_BASE=${array2[0]}
        echo $SAMPLE_ID_BASE	

	# STR Calling
	## EHDN
	echo "Starting EHDN profiling"
	$EHDN8 profile \
		--reads $SAMPLE_BAM \
		--reference $GENOME_FASTA \
		--output-prefix ${SAMPLE_ID_BASE}_Het_EHDN8 

	$EHDN62 --bam $SAMPLE_BAM \
		--reference $GENOME_FASTA \
		--output ${SAMPLE_ID_BASE}_Het_EHDN62.json


done	
	

