# Create pbs scripts to run EHdn on Illumina Polaris cohort

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

PBSDIR=$BASEDIR/scripts/pbs/EHdn_0.8.6_max-irr-mapq60

TEMPLATE=$BASEDIR/scripts/template_EHdn_0.8.6_max-irr-mapq60.sh

# Diversity cohort
mkdir -p $PBSDIR/diversity
echo "cd $PBSDIR/diversity" > $PBSDIR/diversity_submit.sh
for BAM in $(ls $BASEDIR/bam_bwamem/diversity/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/diversity/g" \
        > $PBSDIR/diversity/EHdn_0.8.6_max-irr-mapq60_$BASE.pbs
    echo "qsub EHdn_0.8.6_max-irr-mapq60_$BASE.pbs" >> $PBSDIR/diversity_submit.sh
done


# Create pbs scripts to run RepeatExpansions cohort
mkdir -p $PBSDIR/repeatexpansions
echo "cd $PBSDIR/repeatexpansions" > $PBSDIR/repeatexpansions_submit.sh
for BAM in $(ls $BASEDIR/bam_bwamem/repeatexpansions/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/repeatexpansions/g" \
        > $PBSDIR/repeatexpansions/EHdn_0.8.6_max-irr-mapq60_$BASE.pbs
    echo "qsub EHdn_0.8.6_max-irr-mapq60_$BASE.pbs" >> $PBSDIR/repeatexpansions_submit.sh
done


