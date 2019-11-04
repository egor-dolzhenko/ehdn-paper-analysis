# Variables
CUTOUT_BAM=NA12878_AR_cutout.bam 
REF_BAM=str_bams/AR_Try2/AR_CAG_X-66765158-66765227_X-66763158-66767227_0_25x.sorted.bam
EXP_BAM=str_bams/AR_Try2/AR_CAG_X-66765158-66765227_X-66763158-66767227_138_25x.sorted.bam
OUTPUT_BAM=NA12878_AR_Added_het_138_try2.bam

# merge
samtools merge -r \
	-@ 3 \
	$OUTPUT_BAM \
	$CUTOUT_BAM \
	$REF_BAM \
	$EXP_BAM 
