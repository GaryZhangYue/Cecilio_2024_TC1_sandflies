#!/bin/bash

ml qiime/2-2021.4
# train the taxonomic classifier using SILVA 138 downloaded from https://docs.qiime2.org/2020.6/data-resources with our primer set [341F: 5-CCTAYGGGRBGCASCAG-3; 806R: 5-GGACTACNNGGGTATCTAAT-3)
# download reference files
wget https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza
# extract primer-selected region from full-length 16s database
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTAYGGGRBGCASCAG \
  --p-r-primer GGACTACNNGGGTATCTAAT \
  --o-reads silva-138-99-seqs_341F-806R_extractedReads.qza \
  --p-n-jobs 16 \
  --verbose

# train classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs_341F-806R_extractedReads.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99_341F-806R_classifier.qza

