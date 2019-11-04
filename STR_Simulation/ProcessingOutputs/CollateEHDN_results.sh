#Collate Locus Level
for LocusFile in $(ls *Locus.txt)
do
	sort -k5Vr $LocusFile \
	| head -n 25 \
	> ${LocusFile}.Top25.tsv
done

#Collate Motif Level
for MotifFile in $(ls *Motif_Fixed.txt)
do
        sort -k2Vr $MotifFile \
        | head -n 25 \
        > ${MotifFile}.Top25.tsv
done

cat *Locus.txt.Top25.tsv > LocusCollatedResults.tsv
cat *Motif_Fixed.txt.Top25.tsv > MotifCollatedResults.tsv
