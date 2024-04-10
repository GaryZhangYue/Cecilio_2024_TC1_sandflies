#!/bin/sh

#SBATCH
#SBATCH --job-name=qm2-diversity
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=norm
#SBATCH --mem=200g

        
#load qiime module
ml qiime

# function to signal failure for specific command
function fail {
    echo "FAIL: $@" >&2
        exit 1  # signal failure
        }

source /data/NCBR/projects/NCBR-367/analysis_GYZ/analysis_from_bcdPrm/PE_bcdPrm.config 

echo "Starting moving picture tutorial analysis"
date

echo "Starting alignment and tree"
date

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences ${REPS} \
--p-n-threads 8 \
--o-alignment ${PREFIX}_aligned-rep-seqs.qza \
--o-masked-alignment ${PREFIX}_masked-aligned-rep-seqs.qza \
--o-tree ${PREFIX}_unrooted-tree.qza \
--o-rooted-tree ${PREFIX}_rooted-tree.qza || fail "phylogeny failed"

echo "Starting diversity analysis"
date

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ${PREFIX}_rooted-tree.qza \
  --i-table ${TABLE} \
  --p-sampling-depth ${DEPTH} \
  --m-metadata-file ${METADATA} \
  --output-dir ${PREFIX}_core-metrics-results || fail "core metrics failed"

echo "Starting alpha diversity significance"
date

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${PREFIX}_core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/faith-pd-group-significance.qzv || fail "faith pd failed"

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${PREFIX}_core-metrics-results/shannon_vector.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/shannon-group-significance.qzv || fail "shannon failed"

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${PREFIX}_core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/observed-feature-group-significance.qzv || fail "obs feature failed"


qiime diversity alpha-group-significance \
  --i-alpha-diversity ${PREFIX}_core-metrics-results/evenness_vector.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/evenness-group-significance.qzv || fail "evenness failed"

echo "Beta diversity"
date
qiime emperor plot \
  --i-pcoa ${PREFIX}_core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/weighted-unifrac-emperor.qzv || fail "weighted UniF pcoa failed"

qiime emperor plot \
  --i-pcoa ${PREFIX}_core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/unweighted-unifrac-emperor.qzv || fail "unweighted UniF pcoa failed"

qiime emperor plot \
  --i-pcoa ${PREFIX}_core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_core-metrics-results/bray-curtis-emperor-days.qzv || fail "bray curtis pcoa failed"

echo "Alpha rarefaction"
date

qiime diversity alpha-rarefaction \
  --i-table ${TABLE} \
  --i-phylogeny ${PREFIX}_rooted-tree.qza \
  --p-max-depth ${DEPTH} \
  --m-metadata-file ${METADATA} \
  --o-visualization ${PREFIX}_alpha-rarefaction.qzv || fail "alpha rarefaction failed"

echo "End of script"
date
