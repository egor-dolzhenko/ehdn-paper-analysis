# Compare performance of ExpansionHunter Denovo using different aligners

## Preliminaries

Create directory to store results `BASEDIR`.

Download all scripts here to `$BASEDIR/scripts`.

Download bam file for at least one sample from the [Illumina Polaris Diversity cohort](https://github.com/Illumina/Polaris/wiki/HiSeqX-Diversity-Cohort).

Local configuration settings must be updated in `config.sh`, including:
- path for analysis base directory (`BASEDIR`)
- location of downloaded bam files for Polaris cohorts
- loading particular software (eg `module load XXXX`) if necessary

**Warning:** scripts may include hardcoded paths for local HPC where analysis was performed.


## Running scripts

Run the following scripts in order. Any PBS scripts created need to be successfully run before moving on to the following step.

`1_setup_software_reference.sh`: download reference genome and installs necessary software tools. This file specifies how many samples from the Diversity cohort are aligned by all the tools, which can be adjusted (default = 5).

`2_create_pbs_convert_fastq.sh`: create PBS scripts to convert Isaac aligned [Diversity cohort data](https://github.com/Illumina/Polaris/wiki/HiSeqX-Diversity-Cohort) bam files to fastq. (**Note:** bwamem alignment is performed `polaris_bwameme/` analysis scripts so may not need to be repeated.)

`3_create_pbs_align.sh`: create PBS scripts to align using different tools.

`4_create_pbs_EHdn_0.8.6.sh`: create PBS scripts to run ExpansionHunter Denovo (v0.8.6) across all aligned samples.

