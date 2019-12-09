# Realign Polaris Diversity and RepeatExpansions cohorts with bwamem and run RE tools

Scripts include hardcoded paths for WEHI HPC

Script `1_setup_software_reference.sh` downloads reference genome and installs software tools

Majority of the remaining scripts generate pbs scripts ready for submission to HPC cluster which then need to be submitted.

Local configuration settings must be updated in `config.sh`, including:
- paths for analysis base directory
- location of downloaded bam files for Polaris cohorts
- loading particular software (eg `module load XXXX`)

`BASEDIR` should also be modified in to match value in `config.sh` inside 
`STRetch_pipeline_STRbed_controlsBuiltin.groovy`

Scripts contained this directory are assumed to be copied to $BASEDIR/scripts
