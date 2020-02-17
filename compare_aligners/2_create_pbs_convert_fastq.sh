# Create pbs scripts to convert to fastq files

source config.sh

# *** Select which BAM files to use for aligner comparison ***
# (script below uses first five samples in Diversity cohort)

mkdir -p $BASEDIR/bam_isaac_compare
cd $BASEDIR/bam_isaac_compare

for BAM in $(ls $DIVERSITY_BAMDIR/*.ba* | head -n10)
do
    ln -s $BAM
done


# Create pbs scripts
TEMPLATE=$BASEDIR/scripts/template_convert_fastq.sh
PBS_DIR=$BASEDIR/scripts/pbs/convert_fastq
PBS_SUBMIT=$BASEDIR/scripts/pbs/convert_fastq_submit.sh

mkdir -p $PBS_DIR
echo "cd $PBS_DIR" > $PBS_SUBMIT

for BAM in $(ls $BASEDIR/bam_isaac_compare/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | sed "s#__BAM__#$BAM#g" > $PBS_DIR/convert_fastq_$BASE.pbs
    echo "qsub convert_fastq_$BASE.pbs" >> $PBS_SUBMIT
done

