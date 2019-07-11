#!/bin/sh
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
printf "Start: %s\n" "$(date)"
#
set -x
#
#
# Samples file required for Wilcoxon rank sum statistic.
# This may be very slow on full file. 24h+?
#	--i_sample $DATADIR/gtex_rnaseq_prep_sample.tsv \
#
$cwd/python/gtex_rnaseq_sabv.py \
	--i $DATADIR/gtex_rnaseq_prep_median.tsv \
	--o $DATADIR/gtex_rnaseq_sabv_no-Wilcoxon.tsv \
	-v
#
printf "Done: %s\n" "$(date)"
#
