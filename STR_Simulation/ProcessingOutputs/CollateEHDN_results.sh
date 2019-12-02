# Variables
OUTPUT_DIR=/scratch/st-wasserww-1/EHDN_VERSION_COMPARISON/
cd $OUTPUT_DIR
#Collate Locus Level
for LocusFile in $(ls *Locus.txt)
do
	sort -k5Vr $LocusFile \
	> ${LocusFile}.sorted.tsv
done

#Collate Motif Level
for MotifFile in $(ls *Motif.txt)
do
        sort -k2Vr $MotifFile \
        > ${MotifFile}.sorted.tsv
done
