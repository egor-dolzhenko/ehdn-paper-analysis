# Align using hisat2

#PBS -l nodes=1:ppn=13
#PBS -l mem=10gb
#PBS -l walltime=24:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/bam_compare_aligners/hisat2

# -----------------------------------------------------------------------------

THREADS=12

REF=$BASEDIR/reference/GRCh37-lite

HISAT2=$BASEDIR/bin/hisat2

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

FQ1=$BASEDIR/fastq/${BASE}_R1.fastq.gz
FQ2=$BASEDIR/fastq/${BASE}_R2.fastq.gz

mkdir -p $OUTDIR

$HISAT2 -p $THREADS -x $REF -1 $FQ1 -2 $FQ2 > $OUTDIR/$BASE.sam

samtools sort --threads $THREADS -o $OUTDIR/$BASE.bam $OUTDIR/$BASE.sam
samtools index $OUTDIR/$BASE.bam
rm $OUTDIR/$BASE.sam

