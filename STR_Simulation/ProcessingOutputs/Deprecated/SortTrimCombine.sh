rm EHdnSummaryFile.tsv
for file in $(ls ./*OutlierAnalysis.txt)
do
	echo $file >> EHdnSummaryFile.tsv
	sort $file -k5Vr | head -n 25 >> EHdnSummaryFile.tsv
done

