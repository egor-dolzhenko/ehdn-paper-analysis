# Align using bowtie2 in local mode

#PBS -l nodes=1:ppn=13
#PBS -l mem=10gb
#PBS -l walltime=48:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/bam_compare_aligners/bowtie2_local

# -----------------------------------------------------------------------------

THREADS=12

REF=$BASEDIR/reference/GRCh37-lite

BOWTIE2=$BASEDIR/bin/bowtie2

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

FQ1=$BASEDIR/fastq/${BASE}_R1.fastq.gz
FQ2=$BASEDIR/fastq/${BASE}_R2.fastq.gz

mkdir -p $OUTDIR

$BOWTIE2 -p $THREADS -x $REF -1 $FQ1 -2 $FQ2 --local > $OUTDIR/$BASE.sam

samtools sort --threads $THREADS -o $OUTDIR/$BASE.bam $OUTDIR/$BASE.sam
samtools index $OUTDIR/$BASE.bam
rm $OUTDIR/$BASE.sam

