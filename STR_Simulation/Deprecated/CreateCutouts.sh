## Working directory variables
STR_ANALYSIS_DIR=/sp2work/RICHMOND/STR_Analysis/
FASTADIR=$STR_ANALYSIS_DIR/STR_Simulation/str_fastas/
WRITE_DIR=$STR_ANALYSIS_DIR/STR_Simulation/str_fastas/
WORKING_DIR=$WRITE_DIR
CUTOUTPYTHON=${STR_ANALYSIS_DIR}/STR_Simulation/CreateCutoutBed.py

CHROM_LENGTHS_GRCH37=${STR_ANALYSIS_DIR}STR_Simulation/ReferenceFiles/GRCh37-lite.chrom.sizes

## Fastas to simulate. 
Files=(${FASTADIR}*fasta)

#Varsim variables
ART=/home/richmonp/project/RICHMOND/SIMULATION/varsim-0.8.1/varsim_run/ART/art_bin_VanillaIceCream/art_illumina
READLENGTH=101
DEPTH_HOMO=50
DEPTH_HET=25
DEPTHS=( $DEPTH_HOMO $DEPTH_HET )
# calculated for NA12878 Platinum Bam
FRAGMENT_LENGTH=320
FRAGMENT_STDEV=82

# Load modules for mapping and converting
module load bwa/0.7.15
module load samtools/1.3.1

## Load Modules for GATK
module load nixpkgs/16.09
module load picard/2.1.1 
module load gatk/3.8

#####################


# Phase II: Simulate Reads From Fasta, map reads, and convert them within a loop

cd $WORKING_DIR

Files=(${FASTADIR}*fasta)
for fasta in ${Files[@]}
do
	echo $fasta
	IFS='/' read -a array <<< $fasta
	SampleFasta=${array[-1]}
	IFS='.' read -a array2 <<< "$SampleFasta"
	SAMPLE_ID_BASE=${array2[0]}
	echo $SAMPLE_ID_BASE
	IFS='_' read -a array3 <<< $SAMPLE_ID_BASE
	GENE=${array3[0]}
	
	if [ ! -f ${GENE}*bed ]; then

	python $CUTOUTPYTHON -S $SAMPLE_ID_BASE -C $CHROM_LENGTHS_GRCH37 -O ${SAMPLE_ID_BASE}.GRCh37.cutout.bed
	fi
done

