#!/bin/bash
###
# DOWNLOAD:
###
# Download RNA-Seq files from GTEx portal
###
# RNA-Seq: Gene TPMs
###
# Annotations: Sample Attributes
###
# Annotations: Subject Phenotypes
###
#
cwd="$(pwd)"
#
DATADIR="${cwd}/data"
#
###
# Generate gene xrefs, from Ensembl/Biomart and HUGO.
# Download [protein-coding genes](ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding _gene.txt) from EBI.
# Interactively download BIOMART ENSG2NCBI (human) mapping from Ensembl.org/biomart, with NCBI and HUGO IDs and HUGO symbols. Select only those with Ensembl Protein Family IDs, for protein-encoding genes.
# Get TCRD targets from IDG.
# gtex_gene_xref.R
${cwd}/sh/Go_gtex_GeneXref.sh

###
# Process RNA-Seq data to expression profiles (Exfiles). WARNING: Requires big memory computer (120+GB).
# 
# READ: GTEx Subjects data, 1-row/subject.
# READ: GTEx Samples data, 1-row/sample.
# READ: GTEx RNAseq expression TPM data, 1-row/gene, 1-col/sample.
# READ: gene IDs file, from GTEx/Ensembl/HGNC, via gtex_gene_xref.R.
# REMOVE: samples with Hardy score >2 (prefer healthier).
# REMOVE: samples with high degree of autolysis (self-digestion).
# MERGE: Samples and subjects, to 1-row/sample.
# RESHAPE: RNAseq data from 1-col/sample, to 3 cols: gene, sample, TPM.
# REMOVE: genes in pseudoautosomal regions (PAR) of chromosome Y.
# AGGREGATE: samples, computing median TPM by gene+tissue.
# AGGREGATE: samples, computing median TPM by gene+tissue+sex.
# OUTPUT: median TPMs, 1-row/gene+tissue+sex:
# OUTPUT: expression profiles, 1-row/gene+sex: exfiles_eps.tsv.gz
# Requires ~3hr.
${cwd}/sh/Go_gtex_prep.sh

###
# Compute pairwise similarity coefficients between all profiles.
# Requires ~20GB RAM, ~10 min.
# exfiles_similarity_ruzicka.R
${cwd}/sh/Go_exfiles_sim.sh

###
# Compute pairwise correlation coefficients between all profiles.
# exfiles_similarity_wcorr.R
# Requires ~6hr.
${cwd}/sh/Go_exfiles_cor.sh

###
# Combine correlations and similarity into one file.
# exfiles_similarity_post.py
${cwd}/sh/Go_exfiles_post.sh

###
# Generate final output files for Shiny app and downloads.
${cwd}/R/exfiles_final_files.R


