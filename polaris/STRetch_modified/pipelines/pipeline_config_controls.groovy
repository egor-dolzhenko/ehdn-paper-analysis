// Bpipe pipeline config file
// Paths are relative to the directory the pipeline is running in, so absolute
// paths are recommended.

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8

// For exome pipeline only ***Edit before running the exome pipeline***
EXOME_TARGET="path/to/exome_target_regions.bed"
// Uncomment the line below to run the STRetch installation test, or specify your own
//EXOME_TARGET="SCA8_region.bed"

// For bam pipeline only ***Edit before running if using CRAM input format***
CRAM_REF="path/to/reference_genome_used_to_create_cram.fasta"

// STRetch installation location
STRETCH="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch"

// Paths to tools used by the pipeline
bpipe="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/bpipe"
python="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/python"
goleft="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/goleft"
bedtools="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/bedtools"
bwa="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/bwa"
samtools="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/samtools"
mosdepth="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/mosdepth"
bazam="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/bazam.jar"
picard="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/tools/bin/picard.jar"

// Path to reference data
refdir="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/reference-data"

// Decoy reference assumed to have matching .genome file in the same directory
REF="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/reference-data/hg19.STRdecoys.sorted.fasta"
STR_BED="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed"
DECOY_BED="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/reference-data/STRdecoys.sorted.bed"
// By default, uses other samples in the same batch as a control
CONTROL=""
// Uncomment the line below to use a set of WGS samples as controls, or specify your own
CONTROL="/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/STRetch/reference-data/hg19.PCRfreeWGS_143_STRetch_controls.tsv"

