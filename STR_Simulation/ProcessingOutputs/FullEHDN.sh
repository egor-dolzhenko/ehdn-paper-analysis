#!/bin/bash
#PBS -A dri-wasserman
#PBS -V
#PBS -N EHDN_pipeline
#PBS -m bea
#PBS -V
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=32:mem=120gb
#PBS -l walltime=240:0:0


# This is a workflow for taking a directory full of EHDN Jsons, and one-by-one merging them with the Polaris Background
# Which is then followed up by performing an outlier analysis
# Phillip Richmond
# September 20th, 2019


# Step 0: Make path to EHDN github repo, and other global variables
EHDN=/project/dri-wasserman/RICHMOND/ExpansionHunterDenovo/
SAMPLE_DIR=/scratch/richmonp/SIM/
WORKING_DIR=/scratch/richmonp/SIM/
POLARIS_DIR=/scratch/richmonp/SIM/polaris-ehdn/
GENOME_FASTA=/project/dri-wasserman/GENOME/GRCh37/GRCh37-lite.fa

SAMPLE_JSON=${SAMPLE_ID}_EHDN.json
SAMPLE_BAM=ERR2304566_AFF2_GCC_X-147582151-147582211_X-147580151-147584211_200_18x.sorted.bam
IFS='.' read -a array2 <<< "${SAMPLE_BAM}"
SAMPLE_ID=${array2[0]}

# Step 1: Call STRs using EHDN
$EHDN/precompiled/linux/ExpansionHunterDenovo profile \
	--reads $SAMPLE_BAM \
	--reference $GENOME_FASTA \
	--output ${SAMPLE_ID}_EHDN.json


# Step 2: Combine counts of Polaris + Simulated sample
# A) Make Manifest
IFS='/' read -a array <<< $SAMPLE_JSON
SJSON=${array[-1]}
MANIFEST=${WORKING_DIR}${SAMPLE_ID}_manifest.tsv
rm $MANIFEST




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
$EHDN/precompiled/linux/ExpansionHunterDenovo merge --output-prefix $SAMPLE_ID \
	--reference $GENOME_FASTA \
	--manifest $MANIFEST 	
# For reference, this creates a file with combined counts called: ${SAMPLE_ID}.multisample_profile.json

# Step 3A: Run Outlier Analysis at the locus level
OUTPUT_RESULTS_LOCUS=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Locus.txt
$EHDN/scripts/outlier.py locus  \
	--manifest $MANIFEST \
	--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
	--output $OUTPUT_RESULTS_LOCUS

# Step 3B: Run Outlier Analysis at the motif level
OUTPUT_RESULTS_MOTIF=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Motif.txt
$EHDN/scripts/outlier.py locus  \
        --manifest $MANIFEST \
        --multisample-profile ${SAMPLE_ID}.multisample_profile.json \
        --output $OUTPUT_RESULTS_MOTIF

exit

	
