#!/bin/bash
#############################################################################
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
OFILE_RUZICKA="$DATADIR/gtex_rnaseq_profiles_Ruzicka.tsv"
#
OFILE_WPEARSON="$DATADIR/gtex_rnaseq_profiles_WPearson.tsv"
#
${cwd}/python/exfiles_similarity_post.py \
	--i_cor $OFILE_WPEARSON \
	--i_sim $OFILE_RUZICKA \
	--min_sim 0.7 \
	--min_cor 0.7 \
	--max_anticor "-.7" \
	--o $DATADIR/exfiles_ggc.tsv
#
###
