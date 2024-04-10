qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path pair-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
