#!/bin/sh
#
set -x
#
SRCDATADIR="/home/data/GTEx/data"
DATADIR="data"
#
cwd=$(pwd)
#
DEMO="TRUE"
#
#anno_subj_file="${SRCDATADIR}/GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
anno_subj_file="${SRCDATADIR}/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
#
#anno_samp_file="${SRCDATADIR}/GTEx_v7_Annotations_SampleAttributesDS.txt"
anno_samp_file="${SRCDATADIR}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
#
#rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
#rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz"
#
N_demo=1000
demofile="$DATADIR/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm_DEMO-${N}.gct.gz"
if [ "$DEMO" ]; then
	if [ ! -e "$demofile" ]; then
		gunzip -c $rnaseqfile |head -${N_demo} |gzip -c >$demofile
	fi
	rnaseqfile=${demofile}
fi
#
suffix="_alltissues"
#
# "exfiles_eps.tsv" named for and read by Exfiles app.R.
#
$cwd/python/gtex_rnaseq_prep_app.py \
	--i_subject ${anno_subj_file} \
	--i_sample ${anno_samp_file} \
	--i_rnaseq ${rnaseqfile} \
	--i_gene $DATADIR/gtex_gene_xref.tsv \
	--o_sample $DATADIR/gtex_rnaseq_prep_sample${suffix}.tsv \
	--o_median $DATADIR/gtex_rnaseq_prep_median${suffix}.tsv \
	--o_tissue $DATADIR/gtex_rnaseq_prep_tissues${suffix}.tsv \
	--o_profiles $DATADIR/exfiles_eps${suffix}.tsv \
	--keep_all_tissues \
	-v
#
###
# exfiles_tissue_order.tsv exported from Google sheet.
###
# gtex_gene_xref.tsv written by gtex_gene_map.R, from Ensembl/Biomart and HUGO files.
###
