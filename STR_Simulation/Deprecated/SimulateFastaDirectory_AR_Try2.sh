#!/bin/bash

#SBATCH --account=def-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Phase I: Set up variables 
NSLOTS=$SLURM_JOB_CPUS_PER_NODE
## Working directory variables
PROJECT_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/STR_Analysis/STR_Simulation/
FASTADIR=/home/richmonp/project/RICHMOND/SIMULATION/STR_SIM/STR_Analysis/STR_Simulation/str_fastas/AR_Try2/
WRITE_DIR=/home/richmonp/project/RICHMOND/SIMULATION/STR_SIM/STR_Analysis/STR_Simulation/str_bams/AR_Try2/
WORKING_DIR=$WRITE_DIR
CUTOUTPYTHON=${PROJECT_DIR}CreateCutoutBed.py


## Genome processing variables
BWA_INDEX=/project/projects/def-wyeth/GENOME/UCSC/ucsc_hg19_bwa_index
GENOME_FASTA=/project/projects/def-wyeth/GENOME/UCSC/ucsc.hg19.fasta
CHROM_LENGTHS=/project/projects/def-wyeth/RICHMOND/SIMULATION/STR_SIM/STR_Analysis/STR_Simulation/chrom_end_hg19.bed
TMPDIR=$WORKING_DIR'picardtmp/'
mkdir $TMPDIR

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

module load java
source /project/projects/def-wyeth/TOOLS/SIMULATION_ENVIRONMENT/bin/activate 

cd $WORKING_DIR
for fasta in $(ls ${FASTADIR}*fasta)
do
	echo $fasta
	IFS='/' read -a array <<< $fasta
	SampleFasta=${array[-1]}
	IFS='.' read -a array2 <<< "$SampleFasta"
	SAMPLE_ID_BASE=${array2[0]}
	echo $SAMPLE_ID_BASE
	
	python $CUTOUTPYTHON -S $SAMPLE_ID_BASE -C $CHROM_LENGTHS -O $SAMPLE_ID_BASE.cutout.bed

	for DEPTH in ${DEPTHS[@]}
	do
		echo $DEPTH
		SAMPLE_ID=${SAMPLE_ID_BASE}_${DEPTH}x
		echo $SAMPLE_ID
		$ART -l $READLENGTH \
			-f $DEPTH \
			-o $SAMPLE_ID. \
			-m $FRAGMENT_LENGTH \
			-s $FRAGMENT_STDEV \
			-i $fasta 
		
		### Phase III: Map the simulated reads, and convert to bam
		FASTQR1=$WRITE_DIR$SAMPLE_ID.1.fq
		FASTQR2=$WRITE_DIR$SAMPLE_ID.2.fq
		
		#Map with BWA
		bwa mem $BWA_INDEX -t $NSLOTS -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID.sam
		samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
			| samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
		samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'
		
	done

done
