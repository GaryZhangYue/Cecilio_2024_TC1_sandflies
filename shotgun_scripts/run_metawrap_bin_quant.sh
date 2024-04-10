#!/bin/bash

#SBATCH
#SBATCH --job-name=binquant
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=norm
#SBATCH --mem=300g
#SBATCH --gres=lscratch:500

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh; conda activate metawrap-env; source ~/.bash_profile;

metawrap quant_bins -b ../bin_refinement/metawrap_70_10_bins -o . -a ../metaspades/scaffolds.fasta /data/NCBR/projects/NCBR-367/analysis_shotgun/analysis1/trimmed_reads/S*.fastq

conda deactivate

