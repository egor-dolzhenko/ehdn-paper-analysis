SUMMARYHEADER=SummaryHeader.tsv
SUMMARYHITS=SummarizedHits.tsv
POLARIS_DIR=$PWD

cp $SUMMARYHEADER $SUMMARYHITS

for ANNOTATED_POLARIS in $(ls $POLARIS_DIR/*OutlierAnalysis_Locus_Annotated.tsv)
do
        #echo $ANNOTATED_POLARIS
        # Split the full filepath, only take last file (ignoring directory)
        IFS='/' read -a array1 <<< $ANNOTATED_POLARIS
        APFILE=${array1[-1]}
        # Split the filename, sampleID is first in the _ split 
        IFS='_' read -a array2 <<< $APFILE
        SAMPLE_ID=${array2[0]}
	echo $SAMPLE_ID
	tail -n+2 $ANNOTATED_POLARIS | sed -e "s/^/$SAMPLE_ID\t/g" - >> $SUMMARYHITS
done
