# Create pbs scripts to realign Illumina Polaris cohorts using bwamem

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

PBSDIR=$BASEDIR/scripts/pbs/bwamem_realign

TEMPLATE=$BASEDIR/scripts/template_bwamem_realign.sh

# Diversity cohort
mkdir -p $PBSDIR/diversity
echo "cd $PBSDIR/diversity" > $PBSDIR/diversity_submit.sh
for BAM in $(ls $BASEDIR/bam_isaac/diversity/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/diversity/g" \
        > $PBSDIR/diversity/bwamem_realign_$BASE.pbs
    echo "qsub bwamem_realign_$BASE.pbs" >> $PBSDIR/diversity_submit.sh
done


# Create pbs scripts to run RepeatExpansions cohort
mkdir -p $PBSDIR/repeatexpansions
echo "cd $PBSDIR/repeatexpansions" > $PBSDIR/repeatexpansions_submit.sh
for BAM in $(ls $BASEDIR/bam_isaac/repeatexpansions/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/repeatexpansions/g" \
        > $PBSDIR/repeatexpansions/bwamem_realign_$BASE.pbs
    echo "qsub bwamem_realign_$BASE.pbs" >> $PBSDIR/repeatexpansions_submit.sh
done
