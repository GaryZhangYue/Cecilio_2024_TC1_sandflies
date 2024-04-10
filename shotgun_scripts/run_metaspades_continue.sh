#!/bin/sh

#SBATCH
#SBATCH --job-name=spades_bmt
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=largemem
#SBATCH --mem=1000g
#SBATCH --gres=lscratch:500

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

F=ALL_1.fastq
R=ALL_2.fastq
O=metaspades
mkdir -p $O

date
ml spades/3.15.5
spades.py --continue -o $O || fail "assembly failed"
date
