#!/usr/bin/bash
# Simulation approach to create synthetic STR expansions
# This script will take in a table of pathogenic STRs, and simulate STR fastas and cutout bed files
# Then it will generate reads and map them to a corresponding reference genome
# Phillip Richmond
# August 19th 2019
# Updated October 28th 2019
# This is a modified version to accommodate for complex STR expansions

# 0. Initialize 
# Assumes you have bedtools, pysam, bwa, and python3
# This should be in the setup instructions
conda activate STR_environment

# Variables which you need to set

# Where is the github repo cloned to
STR_ANALYSIS_DIR=/mnt/causes-vnx1/PERSONAL/RICHMOND/STR_Analysis/

# Some variables for script names
SIMFASTA=$STR_ANALYSIS_DIR/STR_Simulation/SimFasta.py
CUTOUTPYTHON=$STR_ANALYSIS_DIR/STR_Simulation/CreateCutoutBed.py

# ART variables, relevant to the simulation you want to perform
ART=art_illumina
READLENGTH=150
DEPTH_HOMO=36
DEPTH_HET=18
DEPTHS=( $DEPTH_HOMO $DEPTH_HET )
# Calculated with samtools stats on the human genome BAM file from Polaris
FRAGMENT_LENGTH=460
FRAGMENT_STDEV=115


##################################################

### GRCh37 ###

# 1. Generate Fasta files of the regions
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
STR_TABLE=$STR_ANALYSIS_DIR/CompareSTRDatabases/PathogenicLoci_GRCh37*table.tsv

python $SIMFASTA \
	-F $GENOME_FASTA \
	-D $STR_ANALYSIS_DIR/STR_Simulation/str_fastas/GRCh37_20191028 \
	-T $STR_TABLE \
	-S 0,10,20,50,100,1000 \
	-W 2000

# 2. Create Cutout BED files for the FASTA files you just generated

FASTADIR=${STR_ANALYSIS_DIR}STR_Simulation/
CHROM_LENGTHS=${STR_ANALYSIS_DIR}STR_Simulation/ReferenceFiles/GRCh37-lite.chrom.sizes
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
       
       if [ ! -f $FASTADIR${GENE}*bed ]; then

       python $CUTOUTPYTHON -S $SAMPLE_ID_BASE -C $CHROM_LENGTHS -O $FASTADIR${SAMPLE_ID_BASE}.GRCh37.cutout.bed
       fi
done

# 3. Simulate the fasta files into reads and map them with BWAmem

## Fastas to simulate. 
Files=(${FASTADIR}*fasta)

cd $WORKING_DIR
for fasta in $(ls ${FASTADIR}*fasta)
do
        echo $fasta
        IFS='/' read -a array <<< $fasta
        SampleFasta=${array[-1]}
        IFS='.' read -a array2 <<< "$SampleFasta"
        SAMPLE_ID_BASE=${array2[0]}
        echo $SAMPLE_ID_BASE

        for DEPTH in ${DEPTHS[@]}
        do
                echo $DEPTH
                SAMPLE_ID=${SAMPLE_ID_BASE}_${DEPTH}x
                echo $SAMPLE_ID
                $ART -l $READLENGTH \
                        -f $DEPTH \
                        -o ${SAMPLE_ID}. \
                        -m $FRAGMENT_LENGTH \
                        -s $FRAGMENT_STDEV \
                        -i $fasta

                ### Phase III: Map the simulated reads, and convert to bam
                FASTQR1=$WRITE_DIR$SAMPLE_ID.1.fq
                FASTQR2=$WRITE_DIR$SAMPLE_ID.2.fq

                #Map with BWA
                bwa mem $BWA_INDEX -t $THREADS -R "@RG\tID:$SAMPLE_ID\tSM:STR_SIMULATED\tPL:illumina" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID.sam
                samtools view -@ $THREADS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
                        | samtools sort - -@ $THREADS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
                samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

        done

done

exit
## NOTE: BELOW IS CURRENTLY DEPRECATED AS OF 2019-10-28
#######################################################

### GRCh38 ###

# 1. Generate Fasta files of the regions
GENOME_FASTA=/sp1work/Genomes/GSC/GRCh38-lite.fa
STR_TABLE=$STR_ANALYSIS_DIR/CompareSTRDatabases/PathogenicLoci_GRCh38*table.tsv

python $SIMFASTA \
        -F $GENOME_FASTA \
        -D $STR_ANALYSIS_DIR/STR_Simulation/str_fastas/GRCh38 \
        -T $STR_TABLE \
        -S 0,10,20,50,100,1000 \
        -W 2000

# 2. Create Cutout BED files for the FASTA files you just generated

FASTADIR=$PWD/str_fastas/GRCh38/
CHROM_LENGTHS=${STR_ANALYSIS_DIR}STR_Simulation/ReferenceFiles/GRCh38-lite.chrom.sizes
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

       if [ ! -f $FASTADIR${GENE}*bed ]; then

       python $CUTOUTPYTHON -S $SAMPLE_ID_BASE -C $CHROM_LENGTHS -O $FASTADIR${SAMPLE_ID_BASE}.GRCh38.cutout.bed
       fi
done

#######################################################

### hg38 ###

# 1. Generate Fasta files of the regions
GENOME_FASTA=/sp1work/Genomes/GSC/hg38-lite.fa
STR_TABLE=$STR_ANALYSIS_DIR/CompareSTRDatabases/PathogenicLoci_hg38*table.tsv

python $SIMFASTA \
        -F $GENOME_FASTA \
        -D $STR_ANALYSIS_DIR/STR_Simulation/str_fastas/hg38 \
        -T $STR_TABLE \
        -S 0,10,20,50,100,1000 \
        -W 2000

# 2. Create Cutout BED files for the FASTA files you just generated

FASTADIR=$PWD/str_fastas/hg38/
CHROM_LENGTHS=${STR_ANALYSIS_DIR}STR_Simulation/ReferenceFiles/hg38.chrom.sizes
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

       if [ ! -f $FASTADIR${GENE}*bed ]; then

       python $CUTOUTPYTHON -S $SAMPLE_ID_BASE -C $CHROM_LENGTHS -O $FASTADIR${SAMPLE_ID_BASE}.hg38.cutout.bed
       fi
done


