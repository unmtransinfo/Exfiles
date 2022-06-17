#!/bin/bash
###
# https://m.ensembl.org/info/data/biomart/biomart_restful.html#biomartperlapi
# https://m.ensembl.org/info/data/biomart/biomart_restful.html#biomartxml

###
# See also:
# http://uswest.ensembl.org/info/data/biomart/biomart_restful.html#wget
# E.g.
# wget -O result.tsv "http://www.ensembl.org/biomart/martservice?query=$(cat biomart_ENSG2NCBI_query.xml)"
###

#
cwd=$(pwd)
DATADIR=${cwd}/data
#
perl $HOME/../app/biomart-perl/scripts/webExample.pl \
	${cwd}/sh/biomart_ENSG2NCBI_query.xml \
	>$DATADIR/biomart_ENSG2NCBI_human.tsv
#

