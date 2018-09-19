#!/bin/bash
#############################################################################
# exfiles_similarity_post.py may be too fat in memory (70GB+).
# So first filter raw files serially.
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
#3GB
OFILE_RUZICKA="$DATADIR/gtex_rnaseq_profiles_Ruzicka.tsv"
#
#22GB
OFILE_WPEARSON="$DATADIR/gtex_rnaseq_profiles_WPearson.tsv"
OFILE_WPEARSON_SLIMMER="$DATADIR/gtex_rnaseq_profiles_WPearson_slimmer.tsv"
#
cat $OFILE_WPEARSON |head -1 >$OFILE_WPEARSON_SLIMMER
#
set -x
#
cat $OFILE_WPEARSON \
	|perl -ne 'while (<>) {$line=$_; s/^.*\t//; if ($_<-.5||$_>.5) {print $line}}' \
	>>$OFILE_WPEARSON_SLIMMER
#
printf "WPearson SLIMMER, FROM: %d TO: %d\n" \
	$(cat $OFILE_WPEARSON |wc -l) \
	$(cat $OFILE_WPEARSON_SLIMMER |wc -l)
#
${cwd}/python/exfiles_similarity_post.py \
	--i_cor $OFILE_WPEARSON_SLIMMER \
	--i_sim $OFILE_RUZICKA \
	--min_sim 0.7 \
	--min_cor 0.7 \
	--max_anticor "-.7" \
	--o $DATADIR/exfiles_ggc.tsv
#
###
