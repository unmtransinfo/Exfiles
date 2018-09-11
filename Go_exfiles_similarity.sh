#!/bin/bash
#############################################################################
# Default correlation min should be -1, i.e. all including anticorrelated.
# Similarity below some threshold typically not useful.
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
IFILE="$DATADIR/gtex_rnaseq_prep_profiles.tsv"
N_IDCOLS="2"
#
#OFILE_TANIMOTO="$DATADIR/gtex_rnaseq_profiles_Tanimoto.tsv"
#TANIMOTO_MIN="0"
#
OFILE_RUZICKA="$DATADIR/gtex_rnaseq_profiles_Ruzicka.tsv"
RUZICKA_MIN="0"
#
#OFILE_PEARSON="$DATADIR/gtex_rnaseq_profiles_Pearson.tsv"
#
OFILE_SPEARMAN="$DATADIR/gtex_rnaseq_profiles_Spearman.tsv"
#
#python/exfiles_similarity.py --i $IFILE \
#	--o_pearson $OFILE_PEARSON
#
#python/exfiles_similarity.py --i $IFILE --n_idcols "$N_IDCOLS" \
#	--o_spearman $OFILE_SPEARMAN
#
#python/exfiles_similarity.py --i $IFILE \
#	--tanimoto_min $TANIMOTO_MIN --o_tanimoto $OFILE_TANIMOTO
#
python/exfiles_similarity.py --i $IFILE --n_idcols "$N_IDCOLS" \
	--ruzicka_min $RUZICKA_MIN --o_ruzicka $OFILE_RUZICKA
#
$(cwd)/R/exfiles_similarity.R
#
