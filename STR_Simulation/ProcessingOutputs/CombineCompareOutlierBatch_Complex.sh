#!/bin/bash
#PBS -A st-wasserww-1
#PBS -V
#PBS -N Complex_EHDN_pipeline_CombineCompareOutlier
#PBS -m bea
#PBS -V
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=4:0:0
# #PBS -J 0-1

# This is a workflow for taking a directory full of EHDN Jsons, and one-by-one merging them with the Polaris Background
# Which is then followed up by performing an outlier analysis
# Phillip Richmond
# September 20th, 2019
# November 1st, 2019

# Step 0: Make path to EHDN github repo, and other global variables
# CHANGE THIS IF YOU ARE USING OUTSIDE OF PHIL
EHDN_DIR=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo-v0.8.0-linux_x86_64/
EHDN080=${EHDN_DIR}/bin/ExpansionHunterDenovo-v0.8.0
EHDN086=/project/st-wasserww-1/TOOLS/newEHDN/ExpansionHunterDenovo
SAMPLE_DIR=/scratch/st-wasserww-1/STR_SIM/Complex_20191122/
WORKING_DIR=/scratch/st-wasserww-1/STR_SIM/Complex_20191122/
POLARIS_DIR_080_MAPQ40=/scratch/st-wasserww-1/POLARIS_JSONS/EHdn_0.8.0/EHdn_0.8.0_max-irr-mapq40_Diversity/
POLARIS_DIR_080_MAPQ60=/scratch/st-wasserww-1/POLARIS_JSONS/EHdn_0.8.0/EHdn_0.8.0_max-irr-mapq60_Diversity/
POLARIS_DIR_086_MAPQ40=/scratch/st-wasserww-1/POLARIS_JSONS/EHdn_0.8.6/EHdn_0.8.6_max-irr-mapq40/diversity/
POLARIS_DIR_086_MAPQ60=/scratch/st-wasserww-1/POLARIS_JSONS/EHdn_0.8.6/EHdn_0.8.6_max-irr-mapq60/diversity/
GENOME_FASTA=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.fa

source activate EHDN

cd $WORKING_DIR
JSON_FILES=( $SAMPLE_DIR/YEATS2_800-192_Het_EHDN086_MAPQ40.str_profile.json )
PBS_ARRAY_INDEX=0
SAMPLE_JSON=${JSON_FILES[$PBS_ARRAY_INDEX]}

# Step 1: Combine counts of Polaris + Simulated sample
# A) Make Manifest
IFS='/' read -a array <<< $SAMPLE_JSON
SJSON=${array[-1]}
SAMPLE_ID=${SJSON::-5}
IFS='_' read -a array2 <<< $SAMPLE_ID
EHDNVERSION=${array2[3]}
MAPQ=${array2[4]}
MANIFEST=${WORKING_DIR}${SAMPLE_ID}_manifest.tsv
rm $MANIFEST

echo $EHDNVERSION
echo $MAPQ
echo $SJSON
# Figure out which POLARIS to use based on the JSON name split
if [ $MAPQ=="MAPQ40.str" ] && [ $EHDNVERSION == "EHDN080" ]
then
	POLARIS_DIR=$POLARIS_DIR_080_MAPQ40
	EHDN=$EHDN080
fi
if [ $MAPQ=="MAPQ60.str" ] && [ $EHDNVERSION == "EHDN080" ] 
then
	POLARIS_DIR=$POLARIS_DIR_080_MAPQ60
	EHDN=$EHDN080
fi

if [ $MAPQ=="MAPQ40.str" ] && [ $EHDNVERSION == "EHDN086" ]
then
	POLARIS_DIR=$POLARIS_DIR_086_MAPQ40
	EHDN=$EHDN086
fi
if [ $MAPQ=="MAPQ60.str" ] && [ $EHDNVERSION == "EHDN086" ]
then
	POLARIS_DIR=$POLARIS_DIR_086_MAPQ60
	EHDN=$EHDN086
fi

echo $POLARIS_DIR

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
OUTPUT_RESULTS_LOCUS=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_EHDN0.8.0_Locus.txt
$PYTHON $EHDN_DIR/scripts/outlier.py locus  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_LOCUS

# And Sort the results
(head -n1 $OUTPUT_RESULTS_LOCUS && tail -n +2 $OUTPUT_RESULTS_LOCUS | sort -k5Vr ) > ${OUTPUT_RESULTS_LOCUS}.sorted.tsv
mv ${OUTPUT_RESULTS_LOCUS}.sorted.tsv $OUTPUT_RESULTS_LOCUS

# Step 2B: Run Outlier Analysis at the motif level
OUTPUT_RESULTS_MOTIF=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_EHDN0.8.0_Motif.txt
$PYTHON $EHDN_DIR/scripts/outlier.py motif  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_MOTIF

# And Sort the results
(head -n1 $OUTPUT_RESULTS_MOTIF && tail -n +2 $OUTPUT_RESULTS_MOTIF | sort -k2Vr ) > ${OUTPUT_RESULTS_MOTIF}.sorted.tsv
mv ${OUTPUT_RESULTS_MOTIF}.sorted.tsv $OUTPUT_RESULTS_MOTIF

