#!/bin/bash

#SBATCH --account=def-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=30G
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Step 0: Load modules for python3 and scipy, make path to EHDN github repo, and other global variables
module load python/3
module load scipy-stack
EHDN=/project/projects/def-wyeth/TOOLS/ExpansionHunterDenovo/
SAMPLE_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/Simulations/
WORKING_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/
POLARIS_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/EHDN/Polaris/polaris-ehdn/

for SAMPLE_JSON in $(ls $SAMPLE_DIR/*json)
do	

	# Step 1: Combine counts of Polaris + Simulated sample
	
	# A) Make Manifest
	IFS='/' read -a array <<< $SAMPLE_JSON
	SJSON=${array[-1]}
	SAMPLE_ID=${SJSON::-5}
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
	
	
	# B) Combine counts using combine_counts.py
	COMBINED_COUNTS=${WORKING_DIR}${SAMPLE_ID}_combinedCounts.json
	python $EHDN/secondary_analysis/combine_counts.py --manifest $MANIFEST \
		--combinedCounts $COMBINED_COUNTS 
	
	# Step 2: Run Outlier Analysis
	OUTPUT_RESULTS=${WORKING_DIR}${SAMPLE_ID}_OutlierAnalysis.txt
	python $EHDN/secondary_analysis/outlier_analysis.py --manifest $MANIFEST \
		--combinedCounts $COMBINED_COUNTS \
		--results $OUTPUT_RESULTS
done	
	

