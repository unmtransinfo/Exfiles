#!/bin/sh

SRCDATADIR="/renci/irods/gtex"
OUTPUT_DIR="${HOME}/expression-profiles/data/output"
mkdir ${OUTPUT_DIR}

./python/gtex_rnaseq_prep_app.py \
	--i_subject $SRCDATADIR/Annotations/GTEx_v7_Annotations_SubjectPhenotypesDS.txt \
	--i_sample $SRCDATADIR/Annotations/GTEx_v7_Annotations_SampleAttributesDS.txt \
	--i_gene ${HOME}/expression-profiles/data/gtex_gene_xref.tsv \
	--i_rnaseq /home/dockeruser/venv/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz \
	--o_sample ${OUTPUT_DIR}/gtex_rnaseq_prep_sample.tsv \
	--o_median ${OUTPUT_DIR}/gtex_rnaseq_prep_median.tsv \
	--o_profiles ${OUTPUT_DIR}/gtex_rnaseq_prep_profiles.tsv \
	--o_tissue ${OUTPUT_DIR}/gtex_rnaseq_prep_tissues.tsv \
	-v
