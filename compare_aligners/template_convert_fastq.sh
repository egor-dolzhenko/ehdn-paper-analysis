# Convert bam to fastq

#PBS -l nodes=1:ppn=2
#PBS -l mem=10gb
#PBS -l walltime=24:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/fastq

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

module load picard-tools/2.17.3

mkdir -p $OUTDIR

java -jar -Xmx10g -XX:ParallelGCThreads=1 /wehisan/general/system/bioinf-software/bioinfsoftware/picard-tools/picard-tools-2.17.3/picard.jar \
    SamToFastq \
    I=$BAM \
    FASTQ=$OUTDIR/${BASE}_R1.fastq.gz \
    SECOND_END_FASTQ=$OUTDIR/${BASE}_R2.fastq.gz

