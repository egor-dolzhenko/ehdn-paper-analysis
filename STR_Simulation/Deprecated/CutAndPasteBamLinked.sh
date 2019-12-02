# Variables
THREADS=8
GENE=AFF2

NATIVE_BAM=/sp2work/RICHMOND/POLARIS/GRCh37/ERR2304564.sorted.bam
STR_ANALYSIS_DIR=/sp2work/RICHMOND/STR_Analysis/
OUTPUT_DIR=$STR_ANALYSIS_DIR/STR_Simulation/FullSimFiles/

CUTOUT_BED=$STR_ANALYSIS_DIR/STR_Simulation/str_bams/
CUTOUT_BAM=${OUTPUT_DIR}/ERR1955529.${GENE}.cutout.bam

SIM_BASE=AFF2_GCC_X-147582151-147582211_X-147580151-147584211
REF_STR_BAM=${STR_ANALYSIS_DIR}/STR_Simulation/str_bams/GRCh37/${SIM_BASE}_0_18x.sorted.bam
EXP_STR_BAM=${STR_ANALYSIS_DIR}/STR_Simulation/str_bams/GRCh37/${SIM_BASE}_200_18x.sorted.bam
OUTPUT_BAM=${STR_ANALYSIS_DIR}/FullSimFiles/ERR1955529.${SIM_BASE}_200_het.bam

# Check files
ls $NATIVE_BAM
ls $REF_STR_BAM
ls $EXP_STR_BAM

# Cut out the BAM given the bed file
samtools view -@ $THREADS -hu -l $CUTOUT_BED \
	$NATIVE_BAM \
	| samtools merge -r \
	-@ $THREADS \
	$OUTPUT_BAM \
	- \
	$REF_BAM \
	$EXP_BAM

exit
# merge
samtools merge -r \
	-@ 3 \
	$OUTPUT_BAM \
	$CUTOUT_BAM \
	$REF_BAM \
	$EXP_BAM 
