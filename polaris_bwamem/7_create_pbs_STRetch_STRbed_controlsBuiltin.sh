# Create pbs scripts to run GangSTR on Illumina Polaris cohort

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

PBSDIR=$BASEDIR/scripts/pbs/STRetch_STRbed_controlsBuiltin

TEMPLATE=$BASEDIR/scripts/template_STRetch_STRbed_controlsBuiltin.sh

# Diversity cohort
mkdir -p $PBSDIR/diversity
echo "cd $PBSDIR/diversity" > $PBSDIR/diversity_submit.sh
for BAM in $(ls $BASEDIR/bam_bwamem/diversity/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/diversity/g" \
        > $PBSDIR/diversity/STRetch_STRbed_controlsBuiltin_$BASE.pbs
    echo "qsub STRetch_STRbed_controlsBuiltin_$BASE.pbs" >> $PBSDIR/diversity_submit.sh
done


# Create pbs scripts to run RepeatExpansions cohort
mkdir -p $PBSDIR/repeatexpansions
echo "cd $PBSDIR/repeatexpansions" > $PBSDIR/repeatexpansions_submit.sh
for BAM in $(ls $BASEDIR/bam_bwamem/repeatexpansions/*.bam)
do
    BASE=$(basename $BAM ".bam")
    
    cat $TEMPLATE | \
        sed "s#__BAM__#$BAM#g" | sed "s/__COHORT__/repeatexpansions/g" \
        > $PBSDIR/repeatexpansions/STRetch_STRbed_controlsBuiltin_$BASE.pbs
    echo "qsub STRetch_STRbed_controlsBuiltin_$BASE.pbs" >> $PBSDIR/repeatexpansions_submit.sh
done

