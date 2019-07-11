#!/bin/sh
#############################################################################
# Samples file required for Wilcoxon rank sum statistic.
# This may be very slow on full file. 24h+?
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
suffix="_alltissues"
#
#	--i_sample $DATADIR/gtex_rnaseq_prep_sample${suffix}.tsv \
#
$cwd/python/gtex_rnaseq_sabv.py \
	--i $DATADIR/gtex_rnaseq_prep_median${suffix}.tsv \
	--o $DATADIR/gtex_rnaseq_sabv${suffix}.tsv \
	-v
#
printf "Done: %s\n" "$(date)"
#
