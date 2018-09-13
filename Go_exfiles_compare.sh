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
###
#OFILE_TANIMOTO="$DATADIR/gtex_rnaseq_profiles_Tanimoto.tsv"
#TANIMOTO_MIN="0"
#
#python/exfiles_similarity.py --i $IFILE \
#	--tanimoto_min $TANIMOTO_MIN --o_tanimoto $OFILE_TANIMOTO
###
#
OFILE_RUZICKA="$DATADIR/gtex_rnaseq_profiles_Ruzicka.tsv"
RUZICKA_MIN="0"
#
#TOO SLOW!
#python/exfiles_similarity.py --i $IFILE --n_idcols "$N_IDCOLS" \
#	--ruzicka_min $RUZICKA_MIN --o_ruzicka $OFILE_RUZICKA
###
N_PROC="8"
rm -f ${OFILE_RUZICKA}_*
python/exfiles_ruzickaMP.py --i $IFILE --n_idcols "$N_IDCOLS" --nproc $N_PROC \
	--ruzicka_min $RUZICKA_MIN --o_ruzicka $OFILE_RUZICKA
head -1 ${OFILE_RUZICKA}_01 >$OFILE_RUZICKA
for f in $(ls ${OFILE_RUZICKA}_*) ; do
	cat $f |sed -e '1d' >>$OFILE_RUZICKA
done
###
#
OFILE_WPEARSON="$DATADIR/gtex_rnaseq_profiles_WPearson.tsv"
#
${cwd}/R/exfiles_similarity_wcorr.R $IFILE $OFILE_WPEARSON
#
###
