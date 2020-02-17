# Align using bwamem

#PBS -l nodes=1:ppn=13
#PBS -l mem=40gb
#PBS -l walltime=30:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/bam_compare_aligners/bwamem

# -----------------------------------------------------------------------------

THREADS=12

REF=$BASEDIR/reference/GRCh37-lite.fa

BWA=$BASEDIR/bin/bwa

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

FQ1=$BASEDIR/fastq/${BASE}_R1.fastq.gz
FQ2=$BASEDIR/fastq/${BASE}_R2.fastq.gz

mkdir -p $OUTDIR

$BWA mem -M -t $THREADS $REF $FQ1 $FQ2 | \
    samtools view -bSu - | \
    samtools sort -o $OUTDIR/$BASE.bam

