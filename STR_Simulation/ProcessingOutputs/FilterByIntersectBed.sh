# The purpose of this script is to filter EHDN_Locus output by a set of exon regions. 
#conda activate STR_environment

# prepare bed file by extending exons +/- 5kb

BEDFILE=/scratch/st-wasserww-1/expressed_brain.sorted.bed
EXTBEDFILE=/scratch/st-wasserww-1/expressed_brain.extended.sorted.bed
GENOMEFILE=/project/st-wasserww-1/GENOME/GRCh37/GRCh37-lite.genome
#bedtools slop -i $BEDFILE -b 5000 -g $GENOMEFILE > $EXTBEDFILE

#exit

for LOCUSFILE in $(ls /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40*Locus.txt)
do
	bedtools intersect \
	  -header \
  	  -a $LOCUSFILE \
	  -b $EXTBEDFILE \
	  -u \
	  > ${LOCUSFILE}.bedfiltered

done



