# Full Comparison Pipeline

# Step 1: 
# Get databases

## STRetch
## hg19.simpleRepeat_period1-6_dedup.sorted.bed comes from the STRetch database package. Look on their website on how to download, or get from here:
## https://figshare.com/s/1a39be9282c90c4860cd
STRETCH_BED=hg19.simpleRepeat_period1-6_dedup.sorted.bed
## GangSTR
#wget https://s3.amazonaws.com/gangstr/hg19/genomewide/hg19_ver13_1.bed.gz
#gunzip hg19_ver13_1.bed.gz
GANGSTR_BED=hg19_ver13_1.bed

## Pathogenic Loci
## Self-curated list of STR loci within reference genome, available for download with this repo
PATHO_BED=PathogenicLoci_LargeSmall_hg19_20191118.bed

# Step 2 
# Compare full databases to each other
## This requires that you install intervene
## conda install -c bioconda intervene
#source activate STR_environment

## Compare with GangSTR13.1
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
        --save-overlaps \
        --figtype pdf \
	--project "FullSTRdbComparison" \
	--title "Full STR Database Comparison" \
	-o ./FullSTRdb_Intervene/ \
	-i  $PATHO_BED $GANGSTR_BED $STRETCH_BED 

## Strict comparison witgh GangSTR13.1
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
        --save-overlaps \
        --figtype pdf \
	--bedtools-options f=0.9,r \
        --project "FullSTRdbComparison_Strict" \
        --title "Full STR Database Comparison, Strict Overlap Threshold" \
        -o ./FullSTRdb_Intervene_Strict/ \
	-i  $PATHO_BED $GANGSTR_BED $STRETCH_BED 



# Step 3
# Filter BED files to only include repeats with motif length 2-6nt
python SCRIPTS/FilterBedForSize.py -T GangSTR -I $GANGSTR_BED -O ${GANGSTR_BED}.filtered.bed --Min 2 --Max 6 
python SCRIPTS/FilterBedForSize.py -T STRetch -I $STRETCH_BED -O ${STRETCH_BED}.filtered.bed --Min 2 --Max 6 


# Step 4:
intervene venn --names=PathogenicLoci,GangSTR,STRetch --type genomic \
	--save-overlaps  \
	--figtype pdf \
        --project "Motifs-2-6nt_STRdbComparison" \
        --title "MotifLen 2-6nt STR Database Comparison" \
	-o ./MotifLen2-6_Intervene \
	-i $PATHO_BED ${GANGSTR_BED}.filtered.bed ${STRETCH_BED}.filtered.bed


