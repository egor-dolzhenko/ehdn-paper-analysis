# Compare STR Databases

> This directory contains analysis for comparing STR databases provided with "genome-wide" tools STRetch and GangSTR.

## Necessary tools:
+ Intervene (https://github.com/asntech/intervene)

## Curated STR expansions
The STR expansions which result in disease have been collated from various sources, and can be found in this file:
[Pathogenic Loci](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/CompareSTRDatabases/PathogenicLoci_GRCh37_20191023_table.tsv)

Gene|PMID/Source|Chrom|Start (0-based)|End (1-based)|Repeat Motif|Pathogenic Lower Bound
-|-|-|-|-|-|-
AFF2|30503517|X|147,582,151|147,582,211|GCC|200
AFF3|24763282|2|100,721,261|100,721,286|CGG|300
AR|29398703|X|66,765,158|66,765,227|CAG|38
ARX|29946432|X|25,031,771|25,031,815|GCG|20
ATN1|29398703|12|7,045,879|7,045,936|CAG|48
ATXN1|29398703|6|16,327,864|16,327,954|CAG|39
ATXN10|29398703|22|46,191,234|46,191,304|ATTCT|280
ATXN2|29398703|12|112,036,753|112,036,822|CAG|33
ATXN3|29398703|14|92,537,353|92,537,386|CAG|55
ATXN7|29398703|3|63,898,360|63,898,390|CAG|37
ATXN8|29398703|13|70,713,515|70,713,560|CTG|80
C11orf80|18160775|11|66,512,293|66,512,320|CGG|500
C9ORF72|29398703|9|27,573,526|27,573,544|GGGGCC|30
CACNA1A|29398703|19|13,318,672|13,318,711|CAG|20
CBL|EH-repo|11|119,076,999|119,077,032|CGG|500
CNBP|29398703|3|128,891,419|128,891,499|CAGG|50
COMP|9887340|19|18,896,839|18,896,859|GAC|6
CSTB|EH-repo|21|45,196,324|45,196,360|CGCGGGGCGGGG|40
DIP2B|EH-repo|12|50,898,784|50,898,805|GGC|500
DMPK|29398703|19|46,273,462|46,273,522|CTG|50
FMR1|29398703|X|146,993,568|146,993,628|CGG|200
FOXL2|12529855|3|138,664,904|138,664,863|GCN|22
FRA10AC1|15203205|10|95,462,280|95,462,305|CGG|200
FXN|29398703|9|71,652,202|71,652,220|GAA|66
GLS|30970188|2|191,745,598|191,745,646|GCA|650
HOXA13|10839976|7|27,239,445|27,239,480|GCN|18
HOXD13|20974300|2|176,957,787|176,957,831|GCN|22
HTT|29398703|4|3,076,603|3,076,660|CAG|36
JPH3|29398703|16|87,637,893|87,637,935|CTG.CAG|41
LOC642361 / NUTM2B-AS1|31332380|10|81,586,141|81,586,160|CGG|90
LRP12|31332380|8|105,601,201|105,601,227|GGC|90
NOP56|30503517|20|2,633,379|2,633,403|GGCCTG|1500
NOTCH2NLC / NBPF19|31332381|1|145,209,323|145,209,344|GGC|90
PAPBN1|9462747|14|23,790,682|23,790,711|GCN|11
PHOX2B|EH-repo|4|41,747,989|41,748,049|GCN|500
PPP2R2B|29398703|5|146,258,290|146,258,320|CAG|43
RUNX2|9182765|6|45,390,488|45,390,538|GCN|27
SOX3|12428212|X|139,586,482|139,586,526|GCN|26
TBP|29398703|6|170,870,994|170,871,105|CAG|43
TCF4|28832669|18|53,253,386|53,253,458|CAG|40
TMEM185A|7874164|X|148,713,420|148,713,437|GCC|300
XYLT1|30554721|16|17,564,764|17,564,779|GCC|500
ZIC2|11285244|13|100,637,703|100,637,747|GCN|25
ZNF713|25196122|7|55,955,294|55,955,332|CGG|85


## Replicating the analysis
To re-run the analysis, simply execute [RunComparison.sh](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/CompareSTRDatabases/RunComparison.sh)

You will need to acquire the STRetch database, or just use the one within this directory.

## Results from the comparison
As you can see, both GangSTR and STRetch together contain all but one of the pathogenic loci within the curated list, and independently they miss some loci. 

![Venn](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/CompareSTRDatabases/FullSTRdb_Intervene/FullSTRdbComparison_venn.png)

## Summarizing findings

> So what did we find in these comparisons?

+ 5 STRs missed by both GangSTR and STRetch

Gene|PMID/Source|Chrom|Start(0-based)|End(1-based)|Repeat_Motif|Pathogenic_Lower_Bound
--|--|--|--|--|--|--
FOXL2|12529855|3|138664904|138664863|GCN|22
NOTCH2NLC/NBPF19|31332381|1|146228800|146228821|GGC|90
HOXA13|10839976|7|27239445|27239480|GCN|18
PAPBN1|9462747|14|23790682|23790711|GCN|11
SOX3|12428212|X|139586482|139586526|GCN|26

+ [Those missed by STRetch](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/CompareSTRDatabases/FullSTRdb_Intervene/sets/110_PathogenicLoci_GangSTR.bed):

Gene|PMID/Source|Chrom|Start(0-based)|End(1-based)|Repeat_Motif|Pathogenic_Lower_Bound
--|--|--|--|--|--|--
COMP|9887340|19|18896839|18896859|GAC|6
CSTB|EH-repo|21|43776443|43776479|CGCGGGGCGGGG|40
DIP2B|EH-repo|12|50505001|50505022|GGC|500
XYLT1|30554721|16|17470907|17470922|GCC|500
LOC642361/NUTM2B-AS1|31332380|10|79826385|79826404|CGG|90


+ What's interesting about this, is there is a paper discussing the discovery of the XYLT1 expansion, detected using STRetch, but the authors were not able to use STRetch "out-of-the-box" to detect the expansion, because the repeat was not defined within the STRetch database.

+ [Those missed by GangSTR](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/CompareSTRDatabases/FullSTRdb_Intervene/sets/101_PathogenicLoci_STRetch.bed):

Gene|PMID/Source|Chrom|Start(0-based)|End(1-based)|Repeat_Motif|Pathogenic_Lower_Bound
--|--|--|--|--|--|--
AFF2|30503517|X|148500631|148500691|GCC|200
AFF3|24763282|2|100721261|100721286|CGG|300
AR|29398703|X|67545316|67545385|CAG|38
ATN1|29398703|12|6936716|6936773|CAG|48
FMR2|29946432|X|148500639|148500685|CCG|200
FXN|29398703|9|69037286|69037304|GAA|66
HOXD13|20974300|2|176957787|176957831|GCN|22
LRP12|31332380|8|104588973|104588999|GGC|90
PHOX2B|EH-repo|4|41745972|41746032|GCN|500
TCF4|28832669|18|55586155|55586227|CAG|40
TMEM185A|7874164|X|148713420|148713437|GCC|300
ZIC2|11285244|13|100637703|100637747|GCN|25

+ It's clear that caution needs to be exerted when using these tools which rely upon predefined repeat loci. 

+ Further, there are repeats expansions which do not exist within the reference sequence. Arguably, they can still be considered repeat expansions because of their ability to expand and contract short tandem repeat motifs, although they are also within a class of non-reference insertions. 

+ Repeats not within the reference include:  

---Non-Reference---
Gene|PMID || Source|Chrom|Start (0-based)|End (1-based)|Repeat Motif|Pathogenic Lower Bound
BEAN1|19878914|16|66,524,302|66,524,369|(TGGAA)exp (TAGAA)n (TAAAA TAGAA)n|500
DAB1|28686858|1|57,832,716|57,832,790|(ATTTT)n (ATTTC)exp (ATTTT)n|31
RAPGEF2|29507423|4|160,263,679|160,263,768|(TTTTA)exp (TTTCA)exp (TTTTA)12|
RFC1|30926972|4|39,350,045|39,350,099|AAGGG|400
SAMD12_a|29507423|8|119379055|119,379,157|(TTTTA)exp (TTTCA)ins|598,458
SAMD12_b|29507423|8|119379055|119,379,157|(TTTTA)exp (TTTCA)ins (TTTTA)exp|2221,225,81
TAF1|29229810|X|70,660,413|70,660,414|CCCTCT |35
TNRC6A|29507423|16|24,624,761|24,624,850|(TTTTA)22 (TTTCA)ins (TTTTA)exp|
YEATS2|31539032|3|183,429,976|183,430,010|(TTTTA)exp (TTTCA)ins|800,192

