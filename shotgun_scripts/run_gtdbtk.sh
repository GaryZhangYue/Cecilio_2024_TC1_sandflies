#!/bin/bash

#SBATCH
#SBATCH --job-name=binClas
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=largemem
#SBATCH --mem=2000g
#SBATCH --gres=lscratch:500

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh
conda activate gtdbtk-2.1.1

mkdir tmp

gtdbtk classify_wf --genome_dir ../bin_refinement/metawrap_70_10_bins/ --out_dir . --cpus 16 --tmpdir ./tmp --full_tree --force --extension fa

conda deactivate

