#!/bin/bash
#############################################################################
# Elapsed ~7m with 18192 eps/sex.
###
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
OFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv"
#
OFILE_WPEARSON="$DATADIR/exfiles_eps_WPearson.tsv"
#
###
printf "Start time: %s\n" "$(date)"
#
${cwd}/python/exfiles_similarity_post.py \
	--i_cor $OFILE_WPEARSON \
	--i_sim $OFILE_RUZICKA \
	--o $DATADIR/exfiles_ggc.tsv
#
printf "Finish time: %s\n" "$(date)"
###
