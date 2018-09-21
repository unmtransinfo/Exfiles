#!/bin/bash
#############################################################################
# 40k profiles can result in 800M results, 80G+ files, so must filter!
# Also low cor/sim values may be uninformative.
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
OFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv.gz"
#
OFILE_WPEARSON="$DATADIR/exfiles_eps_WPearson.tsv.gz"
#
###
# This filtering should be done already in the Ruzicka & WPearson code.
#OFILE_WPEARSON_SLIMMER="$DATADIR/gtex_rnaseq_profiles_WPearson_slimmer.tsv"
#cat $OFILE_WPEARSON |head -1 >$OFILE_WPEARSON_SLIMMER
#cat $OFILE_WPEARSON \
#	|perl -ne 'while (<>) {$line=$_; s/^.*\t//; if ($_<-.5||$_>.5) {print $line}}' \
#	>>$OFILE_WPEARSON_SLIMMER
#printf "WPearson SLIMMER, FROM: %d TO: %d\n" $(cat $OFILE_WPEARSON |wc -l) $(cat $OFILE_WPEARSON_SLIMMER |wc -l)
###
printf "Start time: %s\n" $(date)
#
${cwd}/python/exfiles_similarity_post.py \
	--i_cor $OFILE_WPEARSON \
	--i_sim $OFILE_RUZICKA \
	--min_sim 0.7 \
	--min_cor 0.7 \
	--max_anticor "-.7" \
	--o_unfiltered $DATADIR/exfiles_ggc_unfiltered.tsv \
	--o $DATADIR/exfiles_ggc.tsv
#
printf "Finish time: %s\n" $(date)
###
