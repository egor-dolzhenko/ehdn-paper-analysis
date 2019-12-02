#!/bin/bash
#PBS -A st-wasserww-1
#PBS -V
#PBS -N SetupSTR_Sim_POLARIS
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=32:mem=100gb
#PBS -l walltime=240:0:0

# Phillip Richmond
# August 19 2019

# This script will get reference genomes, index them, and then map a fastq against each of them

###############################################################################
# Step 1 - General Set Up #
###########################

# Some important paths:
# Write DIR for processing
WRITE_DIR=/scratch/st-wasserww-1/

# Where is your genome located
GENOME_DIR=/project/st-wasserww-1/GENOME/

# Where is the github repo cloned to
STR_ANALYSIS_DIR=/project/st-wasserww-1/RICHMOND/STR_Analysis/
FULL_SIM_DIR=$STR_ANALYSIS_DIR/STR_Simulation/FullSimFiles/

# Load necessary tools (should work if you have the conda environment)
source activate STR_environment
# If you haven't made the environment yet, then run:
# conda env create -f $STR_ANALYSIS_DIR/STR_environment.yml

# How many threads will you give this process?
THREADS=30



###############################################################################
# Step 3 - Get Human fastqs, process them to bam #
##################################################

# Get Human Fastq Files
## Polaris project (https://github.com/Illumina/Polaris) HG03522 - AFR - Family:NG109 (Proband,F)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR230/006/ERR2304566/ERR2304566_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR230/006/ERR2304566/ERR2304566_2.fastq.gz

# Process them to BAM
# NOTE: This may be best to do external to this script, since it requires more resources

#### GRCh37 #### 

BWA_INDEX=$GENOME_DIR/GRCh37/GRCh37-lite.fa
GENOME_FASTA=$GENOME_DIR/GRCh37/GRCh37-lite.fa
WORKING_DIR=$WRITE_DIR
FASTQ_DIR=/project/st-wasserww-1/POLARIS/
SAMPLE_ID=ERR2304566

FASTQR1=$FASTQ_DIR${SAMPLE_ID}_1.fastq.gz
FASTQR2=$FASTQ_DIR${SAMPLE_ID}_2.fastq.gz


echo "Working with these read files:"
ls $FASTQR1
ls $FASTQR2

# Map reads
bwa mem $BWA_INDEX \
	-t $THREADS \
	-R "@RG\tID:STR_SIM\tSM:STR_SIM\tPL:illumina" \
	-M \
	$FASTQR1 \
	$FASTQR2 \
	> $WORKING_DIR$SAMPLE_ID.sam


# convert to binary and index
samtools view -@ $THREADS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
	| samtools sort - -@ $THREADS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

#rm $WORKING_DIR$SAMPLE_ID'.sam'

#### GRCh38 ####
BWA_INDEX=$GENOME_DIR/GRCh38/GRCh38-lite.fa
GENOME_FASTA=$GENOME_DIR/GRCh38/GRCh38-lite.fa

FASTQR1=$FASTQ_DIR${SAMPLE_ID}_1.fastq.gz
FASTQR2=$FASTQ_DIR${SAMPLE_ID}_2.fastq.gz


echo "Working with these read files:"
ls $FASTQR1
ls $FASTQR2

# Map reads
bwa mem $BWA_INDEX \
	-t $THREADS \
	-R "@RG\tID:STR_SIM\tSM:STR_SIM\tPL:illumina" \
	-M \
	$FASTQR1 \
	$FASTQR2 \
	> $WORKING_DIR$SAMPLE_ID.sam


# convert to binary and index
samtools view -@ $THREADS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
	| samtools sort - -@ $THREADS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

#rm $WORKING_DIR$SAMPLE_ID'.sam'

