### Config file

# Load required software
module load samtools/1.9
module load java/1.8.0_151

# Base directory where analysis scripts and data are to be stored
BASEDIR=/...path/to/.../ehdn-paper-analysis/compare_aligners

# Directories containing the downloaded Polaris bam files(aligned
# using Isaac) for Diversity cohort
# https://github.com/Illumina/Polaris/wiki/HiSeqX-Diversity-Cohort
DIVERSITY_BAMDIR=/...path/to/.../polaris/bam_isaac/diversity

