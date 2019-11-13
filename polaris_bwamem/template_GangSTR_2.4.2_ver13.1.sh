# Run GangSTR on bwamem realigned Illumina Polaris samples

#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -l walltime=20:00:00

#PBS -m ae
#PBS -M bennett.ma@wehi.edu.au

BAM=__BAM__

BASEDIR=/wehisan/bioinf/lab_bahlo/users/bennett.ma/projects/EHdn/align_bwamem_polaris

OUTDIR=$BASEDIR/REtool_results/GangSTR_2.4.2_ver13.1/__COHORT__

# ---------------------------------------------------------------------------------------

REGIONS=/wehisan/bioinf/lab_bahlo/software/apps/gangSTR/reference_files/hs37_ver13_1.bed

REF=$BASEDIR/reference/GRCh37-lite.fa

module load gcc/6.3.0
GANGSTR=/stornext/HPCScratch/lab_bahlo/Illumina_Polaris/GangSTR/GangSTR-2.4.2/bin/GangSTR

# ---------------------------------------------------------------------------------------

BASE=$(basename $BAM ".bam")

mkdir -p $OUTDIR

# Run GangSTR
$GANGSTR --bam $BAM --ref $REF --regions $REGIONS --out $OUTDIR/$BASE
gzip $OUTDIR/$BASE.vcf

