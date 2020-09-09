#!/bin/bash
#############################################################################
# Elapsed ~7m with 18192 eps/sex.
###
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
IFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv"
IFILE_WPEARSON="$DATADIR/exfiles_eps_WPearson.tsv"
OFILE="$DATADIR/exfiles_ggc.tsv"
###
printf "Start time: %s\n" "$(date)"
#
${cwd}/python/exfiles_similarity_post.py \
	--i_cor $IFILE_WPEARSON \
	--i_sim $IFILE_RUZICKA \
	--o $OFILE
#
printf "Finish time: %s\n" "$(date)"
###
