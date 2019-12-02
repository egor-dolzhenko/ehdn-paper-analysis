# Step 0: Load modules for python3 and scipy, make path to EHDN github repo
module load python/3
module load scipy-stack
EHDN=/project/projects/def-wyeth/TOOLS/ExpansionHunterDenovo/
# Step 1: Combine counts of Polaris + Simulated sample

# A) Make Manifest
WORKING_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/
SAMPLE_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/Simulations/
SAMPLE_ID=pathogenic_str_GLS
SAMPLE_JSON=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/STR_Analysis/STR_Simulation/NA12878_AR_Added_het_138_try2_EHDN.json
POLARIS_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/Polaris/polaris-ehdn/
MANIFEST=${WORKING_DIR}${SAMPLE_ID}_EHDN_manifest.tsv
rm $MANIFEST

echo "$SAMPLE_ID	case	${SAMPLE_JSON}" >> $MANIFEST
for POLARIS_JSON in $(ls $POLARIS_DIR/*json)
do
	# echo $POLARIS_JSON
	# Split the full filepath, only take last file (ignoring directory)
	IFS='/' read -a array <<< $POLARIS_JSON
	PJSON=${array[-1]}
	# Split the filename, removing JSON, could also substring but whatever
	IFS='.' read -a array2 <<< $PJSON
	PJSON_SAMPLE_ID=${array2[0]}
	# echo $PJSON_SAMPLE_ID
	echo "$PJSON_SAMPLE_ID	control	$POLARIS_JSON" >> $MANIFEST
	
done


# B) Combine counts using combine_counts.py
COMBINED_COUNTS=${WORKING_DIR}${SAMPLE_ID}_combinedCounts.json
python $EHDN/secondary_analysis/combine_counts.py --manifest $MANIFEST \
	--combinedCounts $COMBINED_COUNTS 

# Step 2: Run Outlier Analysis
OUTPUT_RESULTS=${WORKING_DIR}${SAMPLE_ID}_EHDN_OutlierAnalysis.txt
python $EHDN/secondary_analysis/outlier_analysis.py --manifest $MANIFEST \
	--combinedCounts $COMBINED_COUNTS \
	--results $OUTPUT_RESULTS

