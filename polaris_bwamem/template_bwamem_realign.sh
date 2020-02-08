# Realign Illumina Polaris samples using bwamem

#PBS -l nodes=1:ppn=13
#PBS -l mem=40gb
#PBS -l walltime=30:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/bam_bwamem/__COHORT__

# -----------------------------------------------------------------------------

THREADS=12

REF=$BASEDIR/reference/GRCh37-lite.fa

BAZAM=$BASEDIR/bin/bazam_1.0.1.jar

BWA=$BASEDIR/bin/bwa-0.7.17/bwa

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

java -Xmx16g -jar \
    $BAZAM -bam $BAM | \
    $BWA mem -M -p -t $THREADS $REF - | \
    samtools view -bSu - | \
    samtools sort -o $OUTDIR/$BASE.bam

