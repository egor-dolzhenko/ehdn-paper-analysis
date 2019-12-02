# Phillip Richmond
# August 19 2019

# This script will get reference genomes, index them, and then map a fastq against each of them

###############################################################################
# Step 1 - General Set Up #
###########################

# Where is the github repo cloned to
STR_ANALYSIS_DIR=/sp2work/RICHMOND/STR_Analysis/
# You can redirect this if you are not simulating within
FULL_SIM_DIR=$STR_ANALYSIS_DIR/STR_Simulation/FullSimFiles/
mkdir -p $FULL_SIM_DIR

# Load necessary tools (should work if you have the conda environment)
conda activate STR_environment
# If you haven't made the environment yet, then run:
# conda env create -f $STR_ANALYSIS_DIR/STR_environment.yml

# How many threads will you give this process?
THREADS=8


# Change to the FullSimFiles directory (note, this directory is ignored by gitignore)
cd $FULL_SIM_DIR

###############################################################################
# Step 2 - Get Reference genomes and index them. #
##################################################

# Get Reference Genome Files (these match the files in str_bams) 
## GRCh37
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai

bwa index GRCh37-lite.fa

## hg38
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip -c Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > GRCh38-lite.fa

bwa index GRCh38-lite.fa


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

BWA_INDEX=$FULL_SIM_DIR/GRCh37-lite.fa
GENOME_FASTA=$FULL_SIM_DIR/GRCh37-lite.fa
WORKING_DIR=$FULL_SIM_DIR
SAMPLE_ID=ERR2304566

FASTQR1=$WORKING_DIR${SAMPLE_ID}_1.fastq.gz
FASTQR2=$WORKING_DIR${SAMPLE_ID}_2.fastq.gz

NSLOTS=8

echo "Working with these read files:"
ls $FASTQR1
ls $FASTQR2

# Map reads
bwa mem $BWA_INDEX \
	-t $NSLOTS \
	-R "@RG\tID:STR_SIM\tSM:STR_SIM\tPL:illumina" \
	-M \
	$FASTQR1 \
	$FASTQR2 \
	> $WORKING_DIR$SAMPLE_ID.sam


# convert to binary and index
samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
	| samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

rm $WORKING_DIR$SAMPLE_ID'.sam'

#### GRCh38 ####

BWA_INDEX=$FULL_SIM_DIR/GRCh38-lite.fa
GENOME_FASTA=$FULL_SIM_DIR/GRCh38-lite.fa
WORKING_DIR=$FULL_SIM_DIR
SAMPLE_ID=ERR2304566

FASTQR1=$WORKING_DIR${SAMPLE_ID}_1.fastq.gz
FASTQR2=$WORKING_DIR${SAMPLE_ID}_2.fastq.gz

NSLOTS=8

echo "Working with these read files:"
ls $FASTQR1
ls $FASTQR2

# Map reads
bwa mem $BWA_INDEX \
	-t $NSLOTS \
	-R "@RG\tID:STR_SIM\tSM:STR_SIM\tPL:illumina" \
	-M \
	$FASTQR1 \
	$FASTQR2 \
	> $WORKING_DIR$SAMPLE_ID.sam


# convert to binary and index
samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
	| samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

rm $WORKING_DIR$SAMPLE_ID'.sam'

