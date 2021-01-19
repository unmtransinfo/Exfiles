#!/bin/bash
###
# https://m.ensembl.org/info/data/biomart/biomart_restful.html#biomartperlapi
# https://m.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml

#
cwd=$(pwd)
DATADIR=${cwd}/data
#
perl $HOME/../app/biomart-perl/scripts/webExample.pl \
	${cwd}/sh/biomart_ENSG2NCBI_query.xml \
	>$DATADIR/biomart_ENSG2NCBI_human.tsv
#
