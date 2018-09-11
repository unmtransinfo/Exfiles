#!/bin/sh

SRCDATADIR="/home/ubuntu/data"
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz"

./python/gtex_rnaseq_prep_app.py \
	--i_subject $SRCDATADIR/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
	--i_sample $SRCDATADIR/GTEx_v7_Annotations_SampleAttributesDS.txt \
	--i_gene $SRCDATADIR/gtex_gene_xref.tsv \
	--i_rnaseq $rnaseqfile \
	--o_sample $SRCDATADIR/gtex_rnaseq_prep_sample.tsv \
	--o_median $SRCDATADIR/gtex_rnaseq_prep_median.tsv \
	--o_profiles $SRCDATADIR/gtex_rnaseq_prep_profiles.tsv \
	--o_tissue $SRCDATADIR/gtex_rnaseq_prep_tissues.tsv \
	-v
