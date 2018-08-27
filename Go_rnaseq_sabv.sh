#!/bin/sh
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
#
#
$cwd/python/gtex_rnaseq_sabv.py \
	--i $DATADIR/gtex_rnaseq_prep_median.tsv \
	--i_gene $HOME/Downloads/biomart_ENSG2NCBI.tsv \
	--o_tau $DATADIR/gtex_rnaseq_sabv_tau.tsv \
	--o_sabv $DATADIR/gtex_rnaseq_sabv.tsv \
	-v
#
#
