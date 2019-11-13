# Run STRetch on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=10
#PBS -l mem=40gb
#PBS -l walltime=20:00:00

#PBS -m ae
#PBS -M bennett.ma@wehi.edu.au

BAM=__BAM__

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

OUTDIR=$BASEDIR/REtool_results/STRetch_STRbed_controlsBuiltin/__COHORT__

# ---------------------------------------------------------------------------------------

STRETCH_DIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris/bin/STRetch

# ---------------------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

CURR_DIR=$PWD

WORK_DIR=$OUTDIR/$BASE
mkdir -p $WORK_DIR
cd $WORK_DIR


# Run STRetch
$STRETCH_DIR/tools/bin/bpipe run -p input_regions=$STRETCH_DIR/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed /wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris/scripts/STRetch_pipeline_STRbed_controlsBuiltin.groovy $BAM





