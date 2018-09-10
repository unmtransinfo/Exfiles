#!/bin/bash
#
### Map RNAseq ENSG IDs to NCBI IDs and HUGO symbols.
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
#
gunzip -c $rnaseqfile \
	|sed -e '1,3d' \
	|awk -F '\t' '{print $1}' \
	|sort -u \
	>$DATADIR/gtex_rnaseq.ensgv
#
printf "Unique VERSIONED_ENSGs: %d\n" $(cat $DATADIR/gtex_rnaseq.ensgv |wc -l)
#
cat $DATADIR/gtex_rnaseq.ensgv \
	|sed -e 's/\..*$//' \
	|sort -u \
	>$DATADIR/gtex_rnaseq.ensg
#
printf "Unique ENSGs: %d\n" $(cat $DATADIR/gtex_rnaseq.ensg |wc -l)
#
#
###
# https://www.genenames.org/cgi-bin/statistics
wget -O - 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt' \
	>$DATADIR/hugo_protein-coding_gene.tsv
###
#API slow. Use HUGO download instead.
#hugo_query.py --get --ftypes "ENSEMBL_GENE_ID" --qfile $DATADIR/gtex_rnaseq.ensg \
#	--o $DATADIR/gtex_rnaseq_ensg_hugo.tsv
###
