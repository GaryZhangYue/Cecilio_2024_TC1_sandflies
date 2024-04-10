#!/bin/bash

ml qiime/2-2021.4
# qiime2 taxonomy classification
source PE_bcdPrm.config

#taxonomy
echo "Starting taxonomic analysis"
date

qiime feature-classifier classify-sklearn \
  --i-classifier ${CLASSI} \
  --i-reads ${REPS} \
  --o-classification ${PREFIX}_taxonomy.qza

qiime metadata tabulate \
  --m-input-file ${PREFIX}_taxonomy.qza \
  --o-visualization ${PREFIX}_taxonomy.qzv

qiime taxa barplot \
	--i-table ${TABLE} \
	--i-taxonomy ${PREFIX}_taxonomy.qza \
	--m-metadata-file ${METADATA} \
	--o-visualization ${PREFIX}_taxonomy_barplots.qzv

