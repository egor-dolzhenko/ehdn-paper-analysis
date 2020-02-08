### Set up software and versions
#
# Use same version of bwa and same reference genome for consistency with other results

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


## bazam version 1.0.1 (https://github.com/ssadedin/bazam/releases/tag/1.0.1)
cd $BASEDIR/bin
curl -L https://github.com/ssadedin/bazam/releases/download/1.0.1/bazam.jar > bazam_1.0.1.jar


## ExpansionHunter Denovo
curl -L https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.8.0/ExpansionHunterDenovo-v0.8.0-linux_x86_64.tar.gz > ExpansionHunterDenovo-v0.8.0-linux_x86_64.tar.gz
tar xzvf ExpansionHunterDenovo-v0.8.0-linux_x86_64.tar.gz
ln -s ExpansionHunterDenovo-v0.8.0-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.0 ExpansionHunterDenovo-v0.8.0
curl -L https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.8.6/ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz > ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz
tar xzvf ExpansionHunterDenovo-v0.8.6-linux_x86_64.tar.gz
ln -s ExpansionHunterDenovo-v0.8.6-linux_x86_64/bin/ExpansionHunterDenovo-v0.8.6 ExpansionHunterDenovo-v0.8.6


## STRetch (use same version for consistency with other results)
git clone https://github.com/Oshlack/STRetch.git
cd $BASEDIR/bin/STRetch
git checkout 5405902
bash install.sh
sed s/chr// reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed > reference-data/GRCh37.simpleRepeat_period1-6_dedup.sorted.bed


