# Index bams

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

PARALLEL_JOBS=40

module load samtools/1.9

# Diversity cohort
parallel -j $PARALLEL_JOBS "samtools index {}" ::: $(ls $BASEDIR/bam_bwamem/diversity/*.bam)

# RepeatExpansions cohort
parallel -j $PARALLEL_JOBS "samtools index {}" ::: $(ls $BASEDIR/bam_bwamem/repeatexpansions/*.bam)


