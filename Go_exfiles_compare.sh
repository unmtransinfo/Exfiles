#!/bin/bash
#############################################################################
#
cwd=$(pwd)
#
DATADIR="$cwd/data"
#
IFILE="$DATADIR/exfiles_eps.tsv"
N_IDCOLS="2"
#
###
#
OFILE_RUZICKA="$DATADIR/exfiles_eps_Ruzicka.tsv"
#
${cwd}/R/exfiles_similarity_ruzicka.R $IFILE $OFILE_RUZICKA
###
#
OFILE_WPEARSON="$DATADIR/exfiles_eps_WPearson.tsv"
#
${cwd}/R/exfiles_similarity_wcorr.R $IFILE $OFILE_WPEARSON
#
###
