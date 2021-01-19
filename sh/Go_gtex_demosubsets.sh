#!/bin/sh
#
SRCDATADIR="/home/data/GTEx/data"
#
cwd=$(pwd)
#
rnaseqfile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
#
#
ofile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo.gct.gz"
#
gunzip -c $rnaseqfile \
	|head -3 \
	|gzip -c >$ofile
#
gunzip -c $rnaseqfile \
	|sed -e '1,3d' \
	|sample_lines.py --i - --p .02 \
	|gzip -c >>$ofile
#
exit
#
#
ofile="$SRCDATADIR/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_demo_10pct.gct.gz"
#
gunzip -c $rnaseqfile \
	|head -3 \
	|gzip -c >$ofile
#
gunzip -c $rnaseqfile \
	|sed -e '1,3d' \
	|sample_lines.py --i - --p .1 \
	|gzip -c >>$ofile
#
#
