# Run EHdn on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=8:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/REtool_results/EHdn_0.8.6_max-irr-mapq60/__COHORT__

# -----------------------------------------------------------------------------

REF=$BASEDIR/reference/GRCh37-lite.fa

EHDN=$BASEDIR/bin/ExpansionHunterDenovo-v0.8.6

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

mkdir -p $OUTDIR

# Run ExpansionHunter Denovo
$EHDN profile --reads $BAM --reference $REF --output-prefix $OUTDIR/$BASE --max-irr-mapq 60


