
source PE_bcdPrm.config

qiime feature-table summarize --i-table ${PREFIX}_dada2.qza --o-visualization ${PREFIX}_dada2_summarize.qzv --m-sample-metadata-file ${METADATA}

qiime feature-table tabulate-seqs \
  --i-data ${PREFIX}_reps.qza \
  --o-visualization ${PREFIX}_reps.qzv

qiime metadata tabulate \
  --m-input-file ${PREFIX}_stats-dada2.qza \
  --o-visualization ${PREFIX}_stats-dada2.qzv
