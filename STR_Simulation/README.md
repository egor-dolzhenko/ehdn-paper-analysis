# STR Simulation

> This subdirectory contains the scripts used for simulation and processing outputs for the ExpansionHunter Denovo paper.


## STR Simulation overview

> Simulation approach used within ExpansionHunter DeNovo paper

This is a multi-step process:
0. Initialize STR simulation environment (BuildSTRSimEnvironment.sh). The simulation framework depends upon:

```
pybedtools
pysam
art
bwa
samtools
htslib
conda
```

1. Create STR fastas, cut-out beds, and then simulate reads and map them to the genome to get STR bams (Create_STR_Simulations\*.sh)
  This is already performed, and repeat expansion fasta and bam files for each str can be found in str_fastas and str_bams respectively. 
2. Get a sample WGS fastq of real data, e.g. from the Polaris Project (GetBackgroundGenome.sh). 
3. Cut out the STR region from the Polaris file, and concatenate with the STR expansion bam file (EmbedAndCallSTRs\*sh). Then call STRs.
4. For processing the output EHdn jsons to outlier prioritized tables, the ProcessingOutputs/ subdirectory uses individual scripts (CombineCompareOutlierBatch\*.sh).
5. Collating these results together uses the script (ParseSTRoutputGeneric.py), which is called with the script (RunParser.sh). 

The final table outputs are used within the supplemental tables S3, S4, S5, S6. Table S7 is curated manually from the OutlierAnalysis files.
