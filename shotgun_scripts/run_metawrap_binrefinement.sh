#!/bin/bash

#SBATCH
#SBATCH --job-name=binref
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=norm
#SBATCH --mem=300g
#SBATCH --gres=lscratch:500

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh; conda activate metawrap-env; source ~/.bash_profile;

metawrap bin_refinement -o . -t 16 -A ../binning/metabat2_bins/ -B ../binning/maxbin2_bins/ -C ../binning/concoct_bins/-c 70 -x 5 || fail "bin refinement failed"

conda deactivate

