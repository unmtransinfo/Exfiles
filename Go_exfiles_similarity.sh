#!/bin/bash
#
cwd=$(pwd)
#
DATADIR=$cwd/data
#
IFILE="$DATADIR/gtex_rnaseq_prep_profiles.tsv"
#
OFILE_TANIMOTO="$DATADIR/gtex_rnaseq_profiles_Tanimoto.tsv"
TANIMOTO_MIN="0.5"
OFILE_ABC="$DATADIR/gtex_rnaseq_profiles_ABC.tsv"
#ABC_MIN="0.0"
OFILE_PEARSON="$DATADIR/gtex_rnaseq_profiles_Pearson.tsv"
PEARSON_MIN="0.5"
OFILE_SPEARMAN="$DATADIR/gtex_rnaseq_profiles_Spearman.tsv"
SPEARMAN_MIN="0.5"
#
python/exfiles_similarity.py --i $IFILE \
	--pearson_min $PEARSON_MIN --o_pearson $OFILE_PEARSON
#
python/exfiles_similarity.py --i $IFILE \
	--spearman_min $SPEARMAN_MIN --o_spearman $OFILE_SPEARMAN
#
python/exfiles_similarity.py --i $IFILE \
	--tanimoto_min $TANIMOTO_MIN --o_tanimoto $OFILE_TANIMOTO
#
python/exfiles_similarity.py --i $IFILE \
	--o_abc $OFILE_ABC
#
