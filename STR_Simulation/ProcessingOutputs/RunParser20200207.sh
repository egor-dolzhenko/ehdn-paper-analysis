#Patho Large
# PathogenicLarge
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLargeEHdn/*EHDN086_MAPQ40_OutlierAnalysis_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLargeEHdn/*EHDN086_MAPQ40_OutlierAnalysis_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLargeSTRetch/*.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20200207_EHDN086_MAPQ40_BedFiltered.tsv

# PathogenicSmall
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmallEHdn/*EHDN086_MAPQ40_OutlierAnalysis_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmallEHdn/*EHDN086_MAPQ40_OutlierAnalysis_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmallSTRetch/*.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20200207_EHDN086_MAPQ40_BedFiltered.tsv

exit



# LOWERBOUND ONLY
# PathogenicLarge
python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191202_EHDN086_MAPQ40_STRetch_LowerBoundOnly.tsv

python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191202_EHDN086_MAPQ60_STRetch_LowerBoundOnly.tsv


# PathogenicSmall
python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191202_EHDN086_MAPQ40_STRetch_LowerBoundOnly.tsv

python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191202_EHDN086_MAPQ60_STRetch_LowerBoundOnly.tsv








# ALL SIZES
# PathogenicLarge
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191202_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191202_EHDN086_MAPQ60_STRetch.tsv






# PathogenicSmall
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191202_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191202_EHDN086_MAPQ60_STRetch.tsv


# NonPathogenic STR
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN086_MAPQ60_STRetch.tsv

# NonPathogenic Large Motif
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN080_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN080_MAPQ60_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_*Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_*Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN086_MAPQ60_STRetch.tsv

