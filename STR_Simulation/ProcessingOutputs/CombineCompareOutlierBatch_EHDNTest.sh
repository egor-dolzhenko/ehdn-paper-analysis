#!/bin/bash
#PBS -A st-wasserww-1
#PBS -V
#PBS -N EHDN_pipeline_CombineCompareOutlier
#PBS -m bea
#PBS -V
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=3:0:0
#PBS -o Sim_STR_output^array_index^.txt
#PBS -e Sim_STR_error^array_index^.txt
#PBS -J 0-11

# This is a workflow for taking a directory full of EHDN Jsons, and one-by-one merging them with the Polaris Background
# Which is then followed up by performing an outlier analysis
# Phillip Richmond
# September 20th, 2019
# November 1st, 2019

# Step 0: Make path to EHDN github repo, and other global variables
# CHANGE THIS IF YOU ARE USING OUTSIDE OF PHIL
EHDN_DIR=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo-v0.8.0-linux_x86_64/
EHDN=${EHDN_DIR}/bin/ExpansionHunterDenovo-v0.8.0
SAMPLE_DIR=/scratch/st-wasserww-1/EHDN_VERSION_COMPARISON/
WORKING_DIR=/scratch/st-wasserww-1/EHDN_VERSION_COMPARISON/
POLARIS_DIR=/scratch/st-wasserww-1/polaris-ehdn/
GENOME_FASTA=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.fa

source activate EHDN

cd $WORKING_DIR
JSON_FILES=( $SAMPLE_DIR/*.json ) 
SAMPLE_JSON=${JSON_FILES[$PBS_ARRAY_INDEX]}



# Step 1: Combine counts of Polaris + Simulated sample
# A) Make Manifest
IFS='/' read -a array <<< $SAMPLE_JSON
SJSON=${array[-1]}
SAMPLE_ID=${SJSON::-5}
MANIFEST=${WORKING_DIR}${SAMPLE_ID}_manifest.tsv
rm $MANIFEST

if [ ! -f ${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Motif.txt ]; then

echo "$SAMPLE_ID	case	${SAMPLE_JSON}" >> $MANIFEST
for POLARIS_JSON in $(ls $POLARIS_DIR/*json)
do
	# echo $POLARIS_JSON
	# Split the full filepath, only take last file (ignoring directory)
	IFS='/' read -a array1 <<< $POLARIS_JSON
	PJSON=${array1[-1]}
	# Split the filename, removing JSON, could also substring but whatever
	IFS='.' read -a array2 <<< $PJSON
	PJSON_SAMPLE_ID=${array2[0]}
	# echo $PJSON_SAMPLE_ID
	echo "$PJSON_SAMPLE_ID	control	$POLARIS_JSON" >> $MANIFEST
done

# B) Combine counts
$EHDN merge --output-prefix $SAMPLE_ID \
--reference $GENOME_FASTA \
--manifest $MANIFEST 	
# For reference, this creates a file with combined counts called: ${SAMPLE_ID}.multisample_profile.json

# Step 2A: Run Outlier Analysis at the locus level
OUTPUT_RESULTS_LOCUS=${WORKING_DIR}${SAMPLE_ID}_oldPolaris_OutlierAnalysis_Locus.txt
$PYTHON $EHDN_DIR/scripts/outlier.py locus  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_LOCUS

# Step 2B: Run Outlier Analysis at the motif level
OUTPUT_RESULTS_MOTIF=${WORKING_DIR}${SAMPLE_ID}_oldPolaris_OutlierAnalysis_Motif.txt
$PYTHON $EHDN_DIR/scripts/outlier.py motif  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_MOTIF
fi

