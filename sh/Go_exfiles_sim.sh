#!/bin/bash
#############################################################################
# Elapsed 17m with 18192 eps/sex.
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
#IFILE="$DATADIR/exfiles_eps.tsv"
IFILE="$DATADIR/exfiles_eps_alltissues.tsv"
#N_IDCOLS="2"
#
###
#
OFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv"
MIN_RUZICKA="0.5"
#
${cwd}/R/exfiles_similarity_ruzicka.R $IFILE $OFILE_RUZICKA $MIN_RUZICKA
###
