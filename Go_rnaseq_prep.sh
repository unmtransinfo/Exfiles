#!/bin/sh
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
# rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz"
#
$cwd/python/gtex_rnaseq_prep_app.py \
	--i_subject $SRCDATADIR/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
	--i_sample $SRCDATADIR/GTEx_v7_Annotations_SampleAttributesDS.txt \
	--i_rnaseq $rnaseqfile \
	--i_gene $DATADIR/biomart_ENSG2NCBI.tsv \
	--o_median $DATADIR/gtex_rnaseq_prep_median.tsv \
	--o_level $DATADIR/gtex_rnaseq_prep_level.tsv \
	-v
#
#
