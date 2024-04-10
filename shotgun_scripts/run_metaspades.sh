#!/bin/sh

#SBATCH
#SBATCH --job-name=spades_bmt
#SBATCH --time=10-00:00:00
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

F=ALL_S_1.fastq
R=ALL_S_2.fastq
O=metaspades
mkdir -p $O

date
ml spades/3.15.5
spades.py --meta -1 $F -2 $R -o $O --threads 16 --memory 960 || fail "assembly failed"
date
