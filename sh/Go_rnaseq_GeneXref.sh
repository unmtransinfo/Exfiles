#!/bin/bash
#############################################################################
### See gtex_gene_xref.R, uses gtex_rnaseq.ensg, which maps
### RNAseq ENSG IDs to NCBI IDs and HUGO symbols.
#############################################################################
#
cwd=$(pwd)
#
SRCDATADIR="$HOME/../data/GTEx/data"
DATADIR="${cwd}/data"
#
set -x
#
# rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
#
# Convert versioned ENSGV to unversioned ENSG:
gunzip -c $rnaseqfile \
	|sed -e '1,3d' \
	|awk -F '\t' '{print $1}' \
	|sort -u \
	>$DATADIR/gtex_rnaseq.ensgv
printf "Unique VERSIONED_ENSGs: %d\n" $(cat $DATADIR/gtex_rnaseq.ensgv |wc -l)
#
cat $DATADIR/gtex_rnaseq.ensgv \
	|sed -e 's/\..*$//' \
	|sort -u \
	>$DATADIR/gtex_rnaseq.ensg
#
printf "Unique ENSGs: %d\n" $(cat $DATADIR/gtex_rnaseq.ensg |wc -l)
###
# gtex_gene_xref.R
# Inputs:
#	ensembl_biomart.tsv (Ensembl.org/biomart, Homo sapiens dataset w/ ENSP IDs.)
#	hugo_protein-coding_gene.tsv (ftp)
#	gtex_rnaseq.ensg
# Output:
#	gtex_gene_xref.tsv (ENSG, NCBI, HGNCID, chr, uniprot, symbol, name)
###
# https://www.genenames.org/cgi-bin/statistics
wget -O - 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt' >$DATADIR/hugo_protein-coding_gene.tsv
###
# Writes gtex_gene_xref.tsv, gtex_gene_xref.uniprot
${cwd}/R/gtex_gene_xref.R
#
###
# IDG:
python3 -m BioClients.idg.Client get_targets \
	--idtype "UNIPROT" \
	--i $DATADIR/gtex_gene_xref.uniprot \
	--o $DATADIR/gtex_gene_idg.tsv
#
