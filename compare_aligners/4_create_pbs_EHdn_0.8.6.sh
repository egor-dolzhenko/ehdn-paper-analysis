# Create pbs scripts to run EHdn for all aligned bams

source config.sh

# Symlink Isaac bam dir
ln -s $BASEDIR/bam_isaac_compare $BASEDIR/bam_compare_aligners/isaac


# Create pbs scripts
TEMPLATE=$BASEDIR/scripts/template_EHdn_0.8.6.sh
PBS_DIR=$BASEDIR/scripts/pbs/EHdn_0.8.6
PBS_SUBMIT=$BASEDIR/scripts/pbs/EHdn_0.8.6_submit.sh

mkdir -p $PBS_DIR
echo "cd $PBS_DIR" > $PBS_SUBMIT

for BAM in $(ls $BASEDIR/bam_compare_aligners/*/*.bam)
do
    BASE=$(basename $BAM ".bam")
    ALIGNER=$(basename $(dirname $BAM))
    
    cat $TEMPLATE | sed "s#__BAM__#$BAM#g" > $PBS_DIR/EHdn_0.8.6_${ALIGNER}_${BASE}.pbs
    echo "qsub $PBS_DIR/EHdn_0.8.6_${ALIGNER}_${BASE}.pbs" >> $PBS_SUBMIT
done


# Adjust EHdn optional parameters for bowtie2 which has maximum MAPQ of 44 
for PBS_SCRIPT in $(ls $PBS_DIR/*bowtie2*)
do
    sed "s#--output-prefix#--min-anchor-mapq 30 --max-irr-mapq 30 --output-prefix#g" $PBS_SCRIPT > ${PBS_SCRIPT}.adjustMAPQ
    mv ${PBS_SCRIPT}.adjustMAPQ $PBS_SCRIPT
done


