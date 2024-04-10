#!/bin/bash

#SBATCH
#SBATCH --job-name=binning
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300g
#SBATCH --gres=lscratch:500

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh; conda activate metawrap-env; source ~/.bash_profile;

metawrap binning  --metabat2 --maxbin2 --concoct -m 2000 -t 16 -a ../metaspades/scaffolds.fasta -o . /data/NCBR/projects/NCBR-367/analysis_shotgun/analysis1/trimmed_reads/S*fastq || fail "binning failed"

conda deactivate

