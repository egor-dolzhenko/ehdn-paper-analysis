# Full Comparison Pipeline

# Step 1: 
# Get databases

## STRetch
## hg19.simpleRepeat_period1-6_dedup.sorted.bed comes from the STRetch database package. Look on their website on how to download, or get from here:
## https://figshare.com/s/1a39be9282c90c4860cd

## GangSTR
#wget https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz
#gunzip hg19_ver13_1.bed.gz

## Pathogenic Loci
## Self-curated list of STR loci within reference genome, available for download with this repo
##  PathogenicLoci_GRCh37_20190708.bed 
PATHOGENIC_LOCI=PathogenicLoci_hg19_20191023_table_RefOnly.bed

# Step 2 
# Compare full databases to each other
## This requires that you install intervene
## conda install -c bioconda intervene

## Compare with GangSTR13.1
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
        --save-overlaps \
        --figtype png \
	--project "FullSTRdbComparison" \
	--title "Full STR Database Comparison" \
	-o ./FullSTRdb_Intervene/ \
	-i  $PATHOGENIC_LOCI hg19_ver13_1.bed  hg19.simpleRepeat_period1-6_dedup.sorted.bed

## Strict comparison witgh GangSTR13.1
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
        --save-overlaps \
        --figtype png \
	--bedtools-options f=0.9,r \
        --project "FullSTRdbComparison_Strict" \
        --title "Full STR Database Comparison, Strict Overlap Threshold" \
        -o ./FullSTRdb_Intervene_Strict/ \
        -i  $PATHOGENIC_LOCI hg19_ver13_1.bed  hg19.simpleRepeat_period1-6_dedup.sorted.bed



# Step 3
# Filter BED files to only include repeats with motif length 2-6nt
python SCRIPTS/FilterBedForSize.py -T GangSTR -I hg19_ver13_1.bed -O hg19_ver13_1_2-6filtered.bed --Min 2 --Max 6 
python SCRIPTS/FilterBedForSize.py -T STRetch -I hg19.simpleRepeat_period1-6_dedup.sorted.bed -O hg19.simpleRepeat_period1-6_dedup.sorted_2-6filtered.bed --Min 2 --Max 6 


# Step 4:
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
	--save-overlaps --bedtools-options f=0.5 \
	--figtype png \
        --project "Motifs-2-6nt_STRdbComparison" \
        --title "MotifLen 2-6nt STR Database Comparison" \
	-o ./MotifLen2-6_Intervene \
	-i $PATHOGENIC_LOCI hg19_ver13_1_2-6filtered.bed hg19.simpleRepeat_period1-6_dedup.sorted_2-6filtered.bed  


