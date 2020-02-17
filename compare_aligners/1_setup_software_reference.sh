### Set up software and versions

source config.sh


## BWA version 0.7.17-r1188 (https://github.com/lh3/bwa/releases/tag/v0.7.17)
mkdir $BASEDIR/bin/
cd $BASEDIR/bin/
curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 > bwa-0.7.17.tar.bz2
tar -jxvf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make

## Reference genome
mkdir $BASEDIR/reference
cd $BASEDIR/reference
wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
# Index
samtools faidx GRCh37-lite.fa
$BASEDIR/bin/bwa-0.7.17/bwa index -a bwtsw GRCh37-lite.fa


## bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
mkdir $BASEDIR/bin/
cd $BASEDIR/bin/
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip/download -O bowtie2-2.3.5.1-linux-x86_64.zip
unzip bowtie2-2.3.5.1-linux-x86_64.zip
ln -s bowtie2-2.3.5.1-linux-x86_64/bowtie2
# Create index
BOWTIE2_DIR=$PWD/bowtie2-2.3.5.1-linux-x86_64
cd $BASEDIR/reference
$BOWTIE2_DIR/bowtie2-build GRCh37-lite.fa GRCh37-lite

## HISAT2 (https://ccb.jhu.edu/software/hisat2/index.shtml)
cd $BASEDIR/bin/
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
ln -s hisat2-2.1.0/hisat2
# Create index
HISAT2_DIR=$PWD/hisat2-2.1.0
cd $BASEDIR/reference
$HISAT2_DIR/hisat2-build GRCh37-lite.fa GRCh37-lite


## ExpansionHunter Denovo
cd $BASEDIR/bin/
curl -L https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.8.6/ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz > ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz
tar xzvf ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz
ln -s ExpansionHunterDenovo-v0.8.6-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.6 ExpansionHunterDenovo-v0.8.6
