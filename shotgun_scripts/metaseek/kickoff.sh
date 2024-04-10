#!/bin/bash
#SBATCH
#SBATCH --job-name=metaseq
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=norm
#SBATCH --mem=16g

ml singularity snakemake python

#python /data/zhangy68/tools/meta-seek/metaseek.py /data/NCBR/projects/NCBR-367/analysis_shotgun/rawdata /data/NCBR/projects/NCBR-367/analysis_shotgun/analysis1 --dry-run
source /data/NCBR/apps/genome-assembly/conda/etc/profile.d/conda.sh
python /data/zhangy68/tools/meta-seek/metaseek.py /data/NCBR/projects/NCBR-367/analysis_shotgun/rawdata /data/NCBR/projects/NCBR-367/analysis_shotgun/analysis1
