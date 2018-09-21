#!/bin/bash
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
IFILE="$DATADIR/exfiles_eps.tsv"
#N_IDCOLS="2"
#
###
#
OFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv"
MIN_RUZICKA="0.5"
#
${cwd}/R/exfiles_similarity_ruzicka.R $IFILE $OFILE_RUZICKA $MIN_RUZICKA
###
