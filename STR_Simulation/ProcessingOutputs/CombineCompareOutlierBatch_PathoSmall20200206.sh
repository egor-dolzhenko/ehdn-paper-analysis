#!/bin/bash
#PBS -A st-wasserww-1
#PBS -V
#PBS -N PathoSmall_EHdnSecondary_CombineCompareOutlier
#PBS -m bea
#PBS -V
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=8:mem=40gb
#PBS -l walltime=2:0:0
#PBS -J 0-140
# This is a workflow for taking a directory full of EHDN Jsons, and one-by-one merging them with the Polaris Background
# Which is then followed up by performing an outlier analysis
# Phillip Richmond
# September 20th, 2019
# November 1st, 2019
# February 3rd, 2020  - This is for testing the null cases

# Step 0: Make path to EHDN github repo, and other global variables
# CHANGE THIS IF YOU ARE USING OUTSIDE OF PHIL
EHDN_DIR=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo/
EHDN=/project/st-wasserww-1/TOOLS/newEHDN/ExpansionHunterDenovo
SAMPLE_DIR=/scratch/st-wasserww-1/STR_SIM/PathoSmallEHdn/
WORKING_DIR=/scratch/st-wasserww-1/STR_SIM/PathoSmallEHdn/
POLARIS_DIR=/scratch/st-wasserww-1/POLARIS_JSONS/EHdn_0.8.6/EHdn_0.8.6_max-irr-mapq40/diversity/
GENOME_FASTA=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.fa

ANNOTATE_EHDN=/project/st-wasserww-1/TOOLS/annotate_EHdn.sh 
ANNOVAR_ANNOTATE_VARIATION=/project/st-wasserww-1/TOOLS/annovar/annotate_variation.pl
ANNOVAR_HUMANDB=/project/st-wasserww-1/TOOLS/annovar/humandb/
ANNOVAR_GENOMEVERSION=hg19

source activate EHDN
cd $WORKING_DIR
JSON_FILES=( $SAMPLE_DIR/*.str_profile.json ) 
SAMPLE_JSON=${JSON_FILES[$PBS_ARRAY_INDEX]}
# used for debugging
#PBS_ARRAY_INDEX=0



# Step 1: Combine counts of Polaris + Simulated sample
# A) Make Manifest
IFS='/' read -a array <<< $SAMPLE_JSON
SJSON=${array[-1]}
IFS='.' read -a array2 <<< $SJSON
SAMPLE_ID=${array2[0]}
IFS='_' read -a array2 <<< $SAMPLE_ID
MANIFEST=${WORKING_DIR}${SAMPLE_ID}_manifest.tsv
rm $MANIFEST


# Used for re-running if failed
if [ -f ${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Locus.txt ]; then
	exit
fi

#echo $SJSON
#echo $POLARIS_DIR
echo $SAMPLE_ID
echo "$SAMPLE_ID	case	${SAMPLE_JSON}" >> $MANIFEST
for POLARIS_JSON in $(ls $POLARIS_DIR/*json)
do
	#echo $POLARIS_JSON
	# Split the full filepath, only take last file (ignoring directory)
	IFS='/' read -a array1 <<< $POLARIS_JSON
	PJSON=${array1[-1]}
	# Split the filename, removing JSON, could also substring but whatever
	IFS='.' read -a array2 <<< $PJSON
	PJSON_SAMPLE_ID=${array2[0]}
	if ! [  $SAMPLE_ID == $PJSON_SAMPLE_ID ];
	then
		echo $PJSON_SAMPLE_ID
		echo $SAMPLE_ID
		echo "$PJSON_SAMPLE_ID	control	$POLARIS_JSON" >> $MANIFEST
	fi
done



# B) Combine counts
$EHDN merge --output-prefix $SAMPLE_ID \
--reference $GENOME_FASTA \
--manifest $MANIFEST 	
# For reference, this creates a file with combined counts called: ${SAMPLE_ID}.multisample_profile.json

# Step 2A: Run Outlier Analysis at the locus level
OUTPUT_RESULTS_LOCUS=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Locus.txt
$PYTHON $EHDN_DIR/scripts/outlier.py locus  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_LOCUS

# And Sort the results
(head -n1 $OUTPUT_RESULTS_LOCUS && tail -n +2 $OUTPUT_RESULTS_LOCUS | sort -k5Vr ) > ${OUTPUT_RESULTS_LOCUS}.sorted.tsv
mv ${OUTPUT_RESULTS_LOCUS}.sorted.tsv $OUTPUT_RESULTS_LOCUS

# Step 2B: Run Outlier Analysis at the motif level
OUTPUT_RESULTS_MOTIF=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis_Motif.txt
$PYTHON $EHDN_DIR/scripts/outlier.py motif  \
--manifest $MANIFEST \
--multisample-profile ${SAMPLE_ID}.multisample_profile.json \
--output $OUTPUT_RESULTS_MOTIF

# And Sort the results
(head -n1 $OUTPUT_RESULTS_MOTIF && tail -n +2 $OUTPUT_RESULTS_MOTIF | sort -k2Vr ) > ${OUTPUT_RESULTS_MOTIF}.sorted.tsv
mv ${OUTPUT_RESULTS_MOTIF}.sorted.tsv $OUTPUT_RESULTS_MOTIF


# Step 3, Annotate with ANNOVAR

sh $ANNOTATE_EHDN \
	--ehdn-results $OUTPUT_RESULTS_LOCUS \
	--ehdn-annotated-results ${OUTPUT_RESULTS_LOCUS::-4}_Annotated.tsv \
	--annovar-annotate-variation $ANNOVAR_ANNOTATE_VARIATION \
	--annovar-humandb $ANNOVAR_HUMANDB \
	--annovar-buildver $ANNOVAR_GENOMEVERSION 

