# Create Diversity custom control database for STRetch
# (https://github.com/Oshlack/STRetch/wiki/Running-STRetch#making-your-own-custom-control-data)

STRETCH_DIR=/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/STRetch/

OUTPUT_DIR=/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/data/STRetch/Diversity/

# Create symlinks from running STRetch on each individual sample to avoid realigning

MERGE_DIR=/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/data/STRetch/Diversity_controls
mkdir -p $MERGE_DIR
cd $MERGE_DIR

ln -s $OUTPUT_DIR/*/*.locus_counts .
ln -s $OUTPUT_DIR/*/*.STR_counts .
ln -s $OUTPUT_DIR/*/*.median_cov .

$STRETCH_DIR/tools/bin/python $STRETCH_DIR/scripts/estimateSTR.py --model $STRETCH_DIR/scripts/STRcov.model.csv --locus_counts *.locus_counts --STR_counts *.STR_counts --median_cov *.median_cov --emit Diversity_controls.tsv

cp Diversity_controls.tsv /stornext/HPCScratch/lab_bahlo/Illumina_Polaris/resources/

