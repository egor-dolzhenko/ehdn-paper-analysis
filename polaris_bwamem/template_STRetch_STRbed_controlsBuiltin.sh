# Run STRetch on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=10
#PBS -l mem=40gb
#PBS -l walltime=20:00:00

BAM=__BAM__

source config.sh

OUTDIR=$BASEDIR/REtool_results/STRetch_STRbed_controlsBuiltin/__COHORT__

# -----------------------------------------------------------------------------

STRETCH_DIR=$BASEDIR/bin/STRetch

# -----------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

CURR_DIR=$PWD

WORK_DIR=$OUTDIR/$BASE
mkdir -p $WORK_DIR
cd $WORK_DIR


# Run STRetch
$STRETCH_DIR/tools/bin/bpipe run -p input_regions=$STRETCH_DIR/reference-data/GRCh37.simpleRepeat_period1-6_dedup.sorted.bed $BASEDIR/scripts/STRetch_pipeline_STRbed_controlsBuiltin.groovy $BAM


