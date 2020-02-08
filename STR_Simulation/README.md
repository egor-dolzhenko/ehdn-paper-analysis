# STR Simulation

> The purpose of this subdirectory is to simulate STRs

## Simulation Pipeline Overview

### STR Simulation embedded within genome 

> Simulation approach used within ExpansionHunter DeNovo paper

This is a multi-step process:
0. Initialize STR simulation environment (BuildSTRSimEnvironment.sh).
1. Create STR fastas, cut-out beds, and then simulate reads and map them to the genome to get STR bams (This is done already for you)
2. Get a sample WGS fastq of real data, e.g. from the Polaris Project (GetBackgroundGenome.sh). 
3. Cut out the STR region from a platinum genome bam file, e.g. Polaris Sample, and concatenate with the STR expansion bam file (). 

### Simulation Diagram

![diagram](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/STR_Simulation/STR_Simulation_Approaches.png)
![diagram]()

### Example for AR repeat expansion, with off-target
#### Repeat Region
![AR_Cutout](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/STR_Simulation/AR_NA12878_Cutout_Plus_Final.png)

#### Off-target hits mapping to TCF4 locus
![OffTarget](https://github.com/Phillip-a-richmond/STR_Analysis/blob/master/STR_Simulation/AR_NA12878_Cutout_Plus_Final_OffTargetHits.png)


## Testing STR Callers

The script CallSTRs_STR_VCFDir.sh works for the option 1 - VarSim-based simulation, on the assumption that the bam files are present.

The script CallSTRs_STR_BAMDir.sh works for option 2 - Embedding-based method. This script will include the cutting and placing of the new region in the bam file. 


### VarSim-based simulation

> This is not used for EHDN paper  

For the VarSim-based simulation, I use the GeneBreaker tool to create STR VCFs, and then simulate whole genome sequences with the STRs added as a VCF within VarSim. 
FullSimulationPipeline_STR_VCFDir.sh - This is a script which can take in the directory of VCFs and run varsim, plus processing, plus STR calling
However, this creates an issue because the STRs which get simulated include all the sequence within the reference genome, which is actually quite different than what we see in a true human dataset. 
For details regarding the VarSim pipeline, find them in the VarSimApproach/ directory.

