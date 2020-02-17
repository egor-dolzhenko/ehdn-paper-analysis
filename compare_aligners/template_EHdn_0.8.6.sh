# Run EHdn on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=4:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/EHdn_0.8.6

# -----------------------------------------------------------------------------

REF=$BASEDIR/reference/GRCh37-lite.fa

EHDN=$BASEDIR/bin/ExpansionHunterDenovo-v0.8.6

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")
ALIGNER=$(basename $(dirname $BAM))

mkdir -p $OUTDIR/$ALIGNER

# Run ExpansionHunter Denovo
$EHDN profile --reads $BAM --reference $REF --output-prefix $OUTDIR/$ALIGNER/$BASE

