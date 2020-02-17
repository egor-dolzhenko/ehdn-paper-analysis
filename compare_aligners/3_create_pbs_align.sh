# Create pbs scripts to align with different tools

source config.sh

JOB_LIST="align_minimap2 align_bowtie2_endtoend align_bowtie2_local align_hisat2"

for JOB in $JOB_LIST
do
    TEMPLATE=$BASEDIR/scripts/template_${JOB}.sh
    PBS_DIR=$BASEDIR/scripts/pbs/$JOB
    PBS_SUBMIT=$BASEDIR/scripts/pbs/${JOB}_submit.sh
    
    mkdir -p $PBS_DIR
    echo "cd $PBS_DIR" > $PBS_SUBMIT
    
    for BAM in $(ls $BASEDIR/bam_isaac_compare/*.bam)
    do
        BASE=$(basename $BAM ".bam")
        
        cat $TEMPLATE | sed "s#__BAM__#$BAM#g" > $PBS_DIR/${JOB}_$BASE.pbs
        echo "qsub ${JOB}_$BASE.pbs" >> $PBS_SUBMIT
    done
done

