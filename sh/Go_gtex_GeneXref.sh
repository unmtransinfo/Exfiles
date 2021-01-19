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
GTEX_VERSION="v8"
GTEX_VERSION_YEAR="2017"
GTEX_VERSION_MONTH="06"
GTEX_VERSION_DAY="05"
rnaseqfile="$SRCDATADIR/GTEx_Analysis_${GTEX_VERSION_YEAR}-${GTEX_VERSION_MONTH}-${GTEX_VERSION_DAY}_${GTEX_VERSION}_RNASeQCv1.1.9_gene_tpm.gct.gz"
#
if [ -e "$rnaseqfile" ]; then
	printf "RNAseq file: \"%s\"\n" "${rnaseqfile}"
else
	printf "Not found: \"%s\"\n" "${rnaseqfile}"
	exit
fi
#
printf "GTEX_VERSION\tGTEX_VERSION_YEAR\tGTEX_VERSION_MONTH\tGTEX_VERSION_DAY\n" \
	>$DATADIR/gtex_info.tsv
printf "${GTEX_VERSION}\t${GTEX_VERSION_YEAR}\t${GTEX_VERSION_MONTH}\t${GTEX_VERSION_DAY}\n" \
	>>$DATADIR/gtex_info.tsv
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
#	biomart_ENSG2NCBI_human.tsv (Ensembl.org/biomart, Homo sapiens dataset w/ ENSP IDs.)
#	hugo_protein-coding_gene.tsv (ftp)
#	gtex_rnaseq.ensg
# Output:
#	gtex_gene_xref.tsv (ENSG, NCBI, HGNCID, chr, uniprot, symbol, name)
###
# https://www.genenames.org/cgi-bin/statistics
wget -O - 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt' >$DATADIR/hugo_protein-coding_gene.tsv
###
#
${cwd}/R/gtex_gene_xref.R \
	 $DATADIR/gtex_rnaseq.ensg \
	 $DATADIR/hugo_protein-coding_gene.tsv \
	 $DATADIR/biomart_ENSG2xrefs_human.tsv \
	 $DATADIR/gtex_gene_xref.tsv
#
###
# IDG:
python3 -m BioClients.idg.tcrd.Client listTargets --o $DATADIR/tcrd_targets.tsv
#
