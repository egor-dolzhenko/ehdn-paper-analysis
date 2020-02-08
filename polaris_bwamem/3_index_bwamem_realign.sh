# Index bams

source config.sh

PARALLEL_JOBS=40

# Diversity cohort
parallel -j $PARALLEL_JOBS "samtools index {}" ::: $(ls $BASEDIR/bam_bwamem/diversity/*.bam)

# RepeatExpansions cohort
parallel -j $PARALLEL_JOBS "samtools index {}" ::: $(ls $BASEDIR/bam_bwamem/repeatexpansions/*.bam)


