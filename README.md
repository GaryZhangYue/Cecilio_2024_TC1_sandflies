# Cecilio 2024 Bioinformatics Analysis Scripts 

This repository contains analysis scripts for the study "Leishmania transmission is disrupted in sandflies colonized by Delftia tsuruhatensis TC1 bacteria". The scripts are organized into two main directories: `16s_scripts` and `shotgun_scripts`.

## Directory Structure

### 1. `16s_scripts`

This directory includes:

- **Bash Scripts:** Scripts to run QIIME2 for generating ASV, taxonomy tables, a phylogenetic tree, and all diversity metrics.
- **R Scripts:** Scripts for downstream analysis.
  - `rarefied_table_analysis_cleanCodes.R`: Analyzes the rarefied table to investigate the relative abundance of genera and diversity metrics.
  - `differentialAbundance_testing_ancombc_cleanCodes.R`: Analyzes the differential abundance of genera.

### 2. `shotgun_scripts`

This directory contains bash scripts used for assembling metagenome-assembled genomes (MAGs) from shotgun metagenomic data. The workflow includes:

- **Preprocessing and Quality Check:** A short Snakemake pipeline (`metaseek`) is used for initial data preprocessing and quality checking.
- **Genome Assembly and Binning:** Scripts for co-assembly, binning, bin refinement, bin quantification, and classification.

