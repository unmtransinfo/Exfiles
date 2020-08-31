#!/bin/bash
#############################################################################
# Elapsed ~6hr with 18192 eps/sex.
###
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
IFILE="$DATADIR/exfiles_eps.tsv"
#N_IDCOLS="2"
#
###
#
OFILE_WPEARSON="$DATADIR/exfiles_eps_WPearson.tsv"
MIN_COR="0.5"
MAX_ANTICOR="-0.5"
#
${cwd}/R/exfiles_similarity_wcorr.R $IFILE $OFILE_WPEARSON $MIN_COR $MAX_ANTICOR
#
###
