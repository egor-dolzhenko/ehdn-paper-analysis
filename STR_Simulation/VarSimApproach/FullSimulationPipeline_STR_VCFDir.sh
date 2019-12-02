#!/bin/bash

#SBATCH --account=def-wyeth

## Mail Options
#SBATCH --mail-user=prichmond@cmmt.ubc.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem=100G
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

## Job Array stuff
#SBATCH --array=0-9%10

# Phase I: Set up variables 
NSLOTS=$SLURM_JOB_CPUS_PER_NODE

## Working directory variables
PROJECT_DIR=/project/projects/def-wyeth/RICHMOND/SIMULATION/
VCFDIR=/home/richmonp/project/RICHMOND/SIMULATION/STR_SIM/str_simulations/NotFinishedToBam/
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

# Simulate the reads using varsim
python ../varsim.py --vc_in_vcf All_20161121.vcf.gz --sv_insert_seq insert_seq.txt \
	--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt \
	--reference  hs37d5.fa --id simu --read_length 150 --vc_num_snp 3000000 --vc_num_ins 100000 \
	--vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 --sv_num_ins 2000 \
	--sv_num_del 2000 --sv_num_dup 200 --sv_num_inv 1000 --sv_percent_novel 0.01 \
	--vc_percent_novel 0.01 --mean_fragment_size 350 --sd_fragment_size 50 \
	--vc_min_length_lim 0 --vc_max_length_lim 49 --sv_min_length_lim 50 \
	--sv_max_length_lim 1000000 --nlanes 1 --total_coverage 50 \
	--simulator_executable ART/art_bin_VanillaIceCream/art_illumina \
	--out_dir $WRITE_DIR${SAMPLE_ID}_out \
	--log_dir $WRITE_DIR${SAMPLE_ID}_log \
	--work_dir $WRITE_DIR${SAMPLE_ID}_work \
	--simulator art \
	--vcfs $VCF

# Rename Fastqs
mv $WRITE_DIR${SAMPLE_ID}_out/lane0.read1.fq.gz $FASTQR1
mv $WRITE_DIR${SAMPLE_ID}_out/lane0.read2.fq.gz $FASTQR2

##########################


## Phase III: Map the simulated reads, and convert to bam
# Load modules for mapping and converting
module load bwa/0.7.15
module load samtools/1.3.1

echo "Primary Analysis Started"
date
cd $WORKING_DIR
#Map with BWA
bwa mem $BWA_INDEX -t $NSLOTS -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID.sam
samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
	| samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

# Load Modules for GATK
module load nixpkgs/16.09
module load picard/2.1.1 
module load gatk/3.8

# Picard MarkDuplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$WORKING_DIR$SAMPLE_ID'.sorted.bam' \
	 REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR \
	 M=$WORKING_DIR$SAMPLE_ID'_DuplicateResults.txt' \
	 O=$WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' 

samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam'

rm $WORKING_DIR$SAMPLE_ID'.sam'
rm $WORKING_DIR$SAMPLE_ID'.sorted.bam'

# GATK Realign
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R $GENOME_FASTA -minReads 5 \
	-I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' \
	-o $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' \
	-nt $NSLOTS 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner \
	-R $GENOME_FASTA -model USE_READS \
	-targetIntervals $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' \
	-I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' \
	-o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'

samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'

# Quick clean up
rm $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' 

echo "Primary Analysis Finished"
date


#################


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
~/project/TOOLS/GangSTR-1.4/bin/GangSTR \
        --bam $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
        --readlength 150 \
        --coverage 50 \
        --ref $GENOME_FASTA \
        --regions $GANGSTR_REGIONS \
        --out ${SAMPLE_ID}_GangSTR

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

# Clean Up
#cp $WRITE_DIR/*vcf $PROJECT_DIR 
#rm $FASTQR1
#rm $FASTQR2



echo "Finished"
date
