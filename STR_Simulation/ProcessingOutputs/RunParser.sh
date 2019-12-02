# LOWER ONLY
# PathogenicLarge
python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191120_EHDN086_MAPQ40_STRetch_LowerBoundOnly.tsv

python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191120_EHDN086_MAPQ60_STRetch_LowerBoundOnly.tsv


# PathogenicSmall
python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN086_MAPQ40_STRetch_LowerBoundOnly.tsv

python ParseSTRoutputGeneric.py \
	-L \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN086_MAPQ60_STRetch_LowerBoundOnly.tsv



# PathogenicLarge
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191120_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Large_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoLarge/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoLarge/*Het.STRs.tsv \
	-O  PathogenicLoci_Large_GRCh37_20191120_EHDN086_MAPQ60_STRetch.tsv


# PathogenicSmall
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN080_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN080_MAPQ60_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Small_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/PathoSmall/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/PathoSmall/*Het.STRs.tsv \
	-O  PathogenicLoci_Small_GRCh37_20191120_EHDN086_MAPQ60_STRetch.tsv
exit





# NonPathogenic STR
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN080_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN080_MAPQ60_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_STR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-S /scratch/st-wasserww-1/STR_SIM/Finished20191115/NonPathoSTR/*Het.STRs.tsv \
	-O  NonPathogenicLoci_STR_GRCh37_20191115_EHDN086_MAPQ60_STRetch.tsv
# NonPathogenic VNTR
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN080_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN080_MAPQ60_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN086_MAPQ40_STRetch.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/NonPathogenicLoci_VNTR_GRCh37_20191115_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Finished20191115/VNTR/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  NonPathogenicLoci_VNTR_GRCh37_20191115_EHDN086_MAPQ60_STRetch.tsv

# Pathogenic Complex
python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Complex_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN080_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  PathogenicLoci_Complex_GRCh37_20191118_EHDN080_MAPQ40.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Complex_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN080_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  PathogenicLoci_Complex_GRCh37_20191118_EHDN080_MAPQ60.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Complex_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN086_MAPQ40.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  PathogenicLoci_Complex_GRCh37_20191118_EHDN086_MAPQ40.tsv

python ParseSTRoutputGeneric.py \
	-T /scratch/st-wasserww-1/STR_Analysis/CompareSTRDatabases/PathogenicLoci_Complex_GRCh37_20191118_table.tsv \
	-EL /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Locus.txt \
	-EM /scratch/st-wasserww-1/STR_SIM/Complex/*EHDN086_MAPQ60.str_profile_OutlierAnalysis_EHDN0.8.0_Motif.txt \
	-O  PathogenicLoci_Complex_GRCh37_20191118_EHDN086_MAPQ60.tsv








