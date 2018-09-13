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
	--i_gene $DATADIR/gtex_gene_xref.tsv \
	--o $DATADIR/exfiles_ggc.tsv
#
###
