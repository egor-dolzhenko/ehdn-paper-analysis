# Create STRetch control file based on Diversity cohort
# (https://github.com/Oshlack/STRetch/wiki/Running-STRetch#making-your-own-custom-control-data)

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

RESULTS_DIR=$BASEDIR/REtool_results/STRetch_STRbed_controlsBuiltin/diversity

MERGEDIR=$BASEDIR/REtool_results/STRetch_diversity_controls

mkdir -p $MERGEDIR
cd $MERGEDIR

ln -s $RESULTS_DIR/*/*.locus_counts .
ln -s $RESULTS_DIR/*/*.STR_counts .
ln -s $RESULTS_DIR/*/*.median_cov .

STRETCH_DIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris/bin/STRetch

$STRETCH_DIR/tools/bin/python $STRETCH_DIR/scripts/estimateSTR.py --model $STRETCH_DIR/scripts/STRcov.model.csv --locus_counts *.locus_counts --STR_counts *.STR_counts --median_cov *.median_cov --emit STRetch_Diversity_control_data.tsv


