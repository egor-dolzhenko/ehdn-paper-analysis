## Working directory variables
WRITE_DIR=/mnt/causes-vnx1/PERSONAL/RICHMOND/STR_Analysis/STR_Simulation/str_jsons/GRCh37/
WORKING_DIR=$WRITE_DIR
mkdir -p $WRITE_DIR

## Genome processing variables
BWA_INDEX=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/GSC/GRCh37-lite.fa
TMPDIR=$WORKING_DIR'picardtmp/'
mkdir $TMPDIR

## VCFs to simulate. These are important since this is what the job array is iterating over
Files=(${VCFDIR}*vcf)
IFS='/' read -a array <<< ${Files[$SLURM_ARRAY_TASK_ID]}
SampleVCF=${array[-1]}
IFS='.' read -a array2 <<< "${SampleVCF}"
SAMPLE_ID=${array2[0]}


