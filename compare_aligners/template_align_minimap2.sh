# Align using minimap2

#PBS -l nodes=1:ppn=13
#PBS -l mem=20gb
#PBS -l walltime=24:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/bam_compare_aligners/minimap2

# -----------------------------------------------------------------------------

THREADS=12

REF=$BASEDIR/reference/GRCh37-lite.fa

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

FQ1=$BASEDIR/fastq/${BASE}_R1.fastq.gz
FQ2=$BASEDIR/fastq/${BASE}_R2.fastq.gz

module load minimap2/2.17

mkdir -p $OUTDIR

minimap2 -ax sr -t $THREADS $REF $FQ1 $FQ2 > $OUTDIR/$BASE.sam

samtools sort --threads $THREADS -o $OUTDIR/$BASE.bam $OUTDIR/$BASE.sam
samtools index $OUTDIR/$BASE.bam
rm $OUTDIR/$BASE.sam

