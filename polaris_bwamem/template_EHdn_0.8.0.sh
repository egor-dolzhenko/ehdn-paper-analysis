# Run EHdn on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=8:00:00

#PBS -m ae
#PBS -M bennett.ma@wehi.edu.au

BAM=__BAM__

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

OUTDIR=$BASEDIR/REtool_results/EHdn_0.8.0/__COHORT__

# ---------------------------------------------------------------------------------------

REF=$BASEDIR/reference/GRCh37-lite.fa

EHDN=/wehisan/bioinf/lab_bahlo/software/apps/ExpansionHunterDenovo/0.8.0/ExpansionHunterDenovo

# ---------------------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

mkdir -p $OUTDIR

# Run ExpansionHunter Denovo
$EHDN profile --reads $BAM --reference $REF --output-prefix $OUTDIR/$BASE


