#!/bin/bash

#SBATCH --account=def-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem-per-cpu=4G
#SBATCH --time=140:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=0-4%5

# Phase I: Set up variables 
NSLOTS=$SLURM_JOB_CPUS_PER_NODE

## Working directory variables
PROJECT_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/
VCFDIR=/home/richmonp/project/RICHMOND/SIMULATION/STR_SIM/str_simulations/RunAllSTR/
WRITE_DIR=/scratch/richmonp/SIMULATION/
WORKING_DIR=$WRITE_DIR

## Genome processing variables
BWA_INDEX=/project/projects/def-wyeth/GENOME/GSC/GRCh37-lite.fa
GENOME_FASTA=/project/projects/def-wyeth/GENOME/GSC/GRCh37-lite.fa
TMPDIR=$WORKING_DIR'picardtmp/'
mkdir $TMPDIR

## VCFs to simulate. These are important since this is what the job array is iterating over
Files=(${VCFDIR}*vcf)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleVCF=${array[-1]}
IFS='.' read -a array2 <<< "${SampleVCF}"
SAMPLE_ID=${array2[0]}
VCF=${Files[$SLURM_ARRAY_TASK_ID]}
echo $VCF
echo $SampleVCF
echo $SAMPLE_ID


#Varsim variables
VARSIM_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/varsim-0.8.1/varsim_run/
FASTQR1=$WRITE_DIR${SAMPLE_ID}_out/$SAMPLE_ID.R1.fq.gz
FASTQR2=$WRITE_DIR${SAMPLE_ID}_out/$SAMPLE_ID.R2.fq.gz

# STR calling variables
GANGSTR_REGIONS=/project/projects/def-wyeth/TOOLS/GangSTR-1.4/data/GRCh37_ver10.sorted.bed
STRETCH_REGIONS=/project/projects/def-wyeth/TOOLS/STRetch/reference-data/grch37.simpleRepeat_period1-6_dedup.sorted.bed 
EH_REGIONS=/project/projects/def-wyeth/TOOLS/ExpansionHunter-v3.0.0-rc2-linux_x86_64/variant_catalog/variant_catalog_grch37.json \


#####################


# Phase II: Simulate Reads From VCF

module load java
source /project/projects/def-wyeth/TOOLS/SIMULATION_ENVIRONMENT/bin/activate 

# Note, because of relative paths here I just jump into the varsim dir for simulation
cd $VARSIM_DIR

### Phase IV: STR Calling
echo "STR Calling Started"
date
cd $WRITE_DIR

# GangSTR
module load intel/2016.4
module load intel/2018.3
module load nlopt
module load gsl
module load htslib
module load gcc/7.3.0

echo $SAMPLE_ID

# STRetch
/project/projects/def-wyeth/TOOLS/STRetch/tools/bpipe-0.9.9.5/bin/bpipe run \
	-p input_regions=$STRETCH_REGIONS \
	/project/projects/def-wyeth/TOOLS/STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy \
	$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'

# ExpansionHunter
/project/projects/def-wyeth/TOOLS/ExpansionHunter-v3.0.0-rc2-linux_x86_64/bin/ExpansionHunter \
	--reads $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
	--reference $GENOME_FASTA \
	--variant-catalog $EH_REGIONS \
	--output-prefix ${SAMPLE_ID}_EH3

# ExpansionHunter Denovo
/project/projects/def-wyeth/TOOLS/ExpansionHunter_DeNovo/ExpansionHunterDenovo-v0.6.2 \
	--bam $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
	--reference  $GENOME_FASTA \
	--output ${SAMPLE_ID}_EHDN.json \
# GangSTR
/project/projects/def-wyeth/TOOLS/GangSTR-1.4/bin/GangSTR \
        --bam $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
        --readlength 150 \
        --coverage 50 \
        --ref $GENOME_FASTA \
        --regions $GANGSTR_REGIONS \
        --out ${SAMPLE_ID}_GangSTR

# Clean Up
#cp $WRITE_DIR/*vcf $PROJECT_DIR 
#rm $FASTQR1
#rm $FASTQR2



echo "Finished"
date
