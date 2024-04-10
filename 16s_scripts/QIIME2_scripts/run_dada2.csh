#!/bin/sh

#SBATCH
#SBATCH --job-name=qiime-dd2
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=norm
#SBATCH --mem=100g

source /data/NCBR/projects/NCBR-367/analysis_GYZ/analysis_from_bcdPrm/PE_bcdPrm.config

ml qiime
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $DEMUX \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --o-representative-sequences ${PREFIX}_reps.qza \
  --o-table ${PREFIX}_dada2.qza \
  --o-denoising-stats ${PREFIX}_stats-dada2.qza \
  --p-n-threads 0 \
  --p-min-fold-parent-over-abundance 10
