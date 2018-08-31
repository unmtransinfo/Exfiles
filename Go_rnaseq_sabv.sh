#!/bin/sh
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
#
$cwd/python/gtex_rnaseq_sabv.py \
	--i $DATADIR/gtex_rnaseq_prep_median.tsv \
	--i_sample $DATADIR/gtex_rnaseq_prep_sample.tsv \
	--o $DATADIR/gtex_rnaseq_sabv.tsv \
	-v
#
#
