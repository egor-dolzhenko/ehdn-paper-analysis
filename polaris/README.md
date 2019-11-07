# Run EHDN, GangSTR and STRetch on Illumina Polaris Diversity and RepeatExpansions cohort

Scripts include hardcoded paths for WEHI HPC

Results generated in `data/` directory copied to ehdn-paper-data on Box

Reference genomes GRCh37 and hg19 downloaded from:
https://support.illumina.com/sequencing/sequencing_software/igenome.html

Directory `STRetch_modified` contains modified / additional files only copied from <STRetch_install_dir>

Software versions:
- ExpansionHunter Denovo 0.8.0
- GangSTR 2.4.2 with database version 13.1
- STRetch commit "5405902" from 24 January (`git checkout 5405902`)

STRetch scripts written to be run in succession to "reuse" BAM files and avoid realigning.
1. pbs scripts run each sample using inbuilt control cohort
2. generate STRetch control cohort file from Diversity cohort
3. rerun using Diversity cohort as control group (**TODO**)
   (is this valid for Diversity cohort itself? would outlier mode without controls be better?)

