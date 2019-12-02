# Phillip Richmond
# Optional Header info for simulating using PBS
#!/bin/bash
#PBS -A st-wasserww-1
#PBS -V
#PBS -N STR_Surgery_Pathogenic_Small
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=32:mem=150gb
#PBS -l walltime=32:0:0
#PBS -J 0-20

# Variables
# CHANGE THIS: Tool paths. If you install tools in specific places not in your path, or use module systems then modify this accordingly
source activate STR_environment
EHDN080=/project/st-wasserww-1/TOOLS/ExpansionHunterDenovo-v0.8.0-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.0 
EHDN086=/project/st-wasserww-1/TOOLS/newEHDN/ExpansionHunterDenovo
GANGSTR=/project/st-wasserww-1/TOOLS/GangSTR-2.4.2/bin/GangSTR
GANGSTR_REGIONS=/project/st-wasserww-1/TOOLS/GangSTR-2.4.2/GRCh37_ver13_1.bed
STRETCH=/project/st-wasserww-1/TOOLS/STRetch/tools/bin/bpipe
STRETCH_DIR=/project/st-wasserww-1/TOOLS/STRetch/
GENOME_FASTA=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.fa
THREADS=30

# CHANGE THIS: Make sure you set this to the POLARIS BAM file you are using for background simulation
NATIVE_BAM=/project/st-wasserww-1/POLARIS/ERR2304597_GRCh37.sorted.bam

# CHANGE THIS: This is the directory where you've cloned STR_Analysis
STR_ANALYSIS_DIR=/scratch/st-wasserww-1/STR_Analysis/

# CHANGE THIS: Where do you want the simulation to occur. Set this to a directory with ample disk space, since each simulated/cut-out-and-replaced BAM is roughly 60G
OUTPUT_DIR=/scratch/st-wasserww-1/STR_SIM/
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Set the CUTOUT, and STR_BAM Directory based on where this repo is cloned
CUTOUT_BED_DIR=$STR_ANALYSIS_DIR/STR_Simulation/str_fastas/GRCh37_20191118/
STR_BAM_DIR=$STR_ANALYSIS_DIR/STR_Simulation/str_bams/GRCh37_20191118/

# This is the STR table you simulated from to create the STR BAMs. I also take it and make it easier to parse below
STR_TABLE=$STR_ANALYSIS_DIR/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv
tail -n+2 $STR_TABLE | sed -e 's/\t/,/g' > ${STR_TABLE}.noheader

# Map this Table into an array, so I can loop over it using PBS
mapfile -t myArray2 < ${STR_TABLE}.noheader

# Next I'll parse out the gene, and pathogenic lower bound.
# This is running in a PBS Job Array, where the index will be the row number from the table 
STRline=${myArray2[$PBS_ARRAY_INDEX]}
echo $STRline
IFS=',' read -ra STRColumns <<< "$STRline"
Gene=${STRColumns[0]}
Size1=${STRColumns[6]}

# I'll also expand these sizes
Size2=$((Size1 + 10))
Size3=$((Size1 + 20))
Size4=$((Size1 + 50))
Size5=$((Size1 + 100))
Size6=$((Size1 + 1000))

# Now I'll loop through the sizes, perform the surgery on the BAM files
# And then call STRs
Sizes=( $Size1 $Size2 $Size3 $Size4 $Size5 $Size6 )
for Size in ${Sizes[@]} 
do
	REF_STR_BAM=${STR_BAM_DIR}/${Gene}*_0_18x.sorted.bam
	EXP_STR_BAM=${STR_BAM_DIR}/${Gene}*_${Size}_18x.sorted.bam
	OUTPUT_BAM=$OUTPUT_DIR/${Gene}_${Size}_Het.bam
	CUTOUT_BED=$CUTOUT_BED_DIR/$Gene*bed
	ls $REF_STR_BAM
	ls $EXP_STR_BAM
	ls $CUTOUT_BED
	echo $OUTPUT_BAM
	
	# BAM Surgery
	## Cut-out, then merge with wild-type and expanded BAMs
	if [ ! -f ${Gene}_${Size}_Het.STRs.tsv ]; then
	echo "Starting BAM Surgery"
	samtools view -@ $THREADS -hu -L $CUTOUT_BED \
		$NATIVE_BAM \
		| samtools merge -r \
		-@ $THREADS \
		$OUTPUT_BAM \
		- \
		$REF_STR_BAM \
		$EXP_STR_BAM
	## Index
	samtools index $OUTPUT_BAM

	# STR Calling
	## EHDN
	echo "Starting EHDN profiling"
	$EHDN080 profile \
		--reads $OUTPUT_BAM \
		--reference $GENOME_FASTA \
		--output-prefix ${Gene}_${Size}_Het_EHDN080_MAPQ40 

	$EHDN080 profile \
		--reads $OUTPUT_BAM \
		--reference $GENOME_FASTA \
		--max-irr-mapq 60 \
		--output-prefix ${Gene}_${Size}_Het_EHDN080_MAPQ60 

	$EHDN086 profile \
		--reads $OUTPUT_BAM \
		--reference $GENOME_FASTA \
		--output-prefix ${Gene}_${Size}_Het_EHDN086_MAPQ40 

	$EHDN086 profile \
		--reads $OUTPUT_BAM \
		--reference $GENOME_FASTA \
		--max-irr-mapq 60 \
		--output-prefix ${Gene}_${Size}_Het_EHDN086_MAPQ60 

	## STRetch
	echo "Starting STRetch"
	$STRETCH run \
		-p input_regions=$STRETCH_DIR/reference-data/GRCh37.simpleRepeat_period1-6_dedup.sorted.bed \
		$STRETCH_DIR/pipelines/STRetch_wgs_bam_pipeline.groovy \
		$OUTPUT_BAM

	## GangSTR
#	echo "Starting GangSTR"
#	$GANGSTR  \
#		--bam $OUTPUT_BAM \
#		--ref $GENOME_FASTA \
#		--regions $GANGSTR_REGIONS \
#		--out ${Gene}_${Size}_Het_GangSTR
#
	## Clean Up
	rm $OUTPUT_BAM
	fi
done

