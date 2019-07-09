#!/bin/bash

cwd=$(pwd)

DATADIR="${cwd}/data"

${cwd}/sh/runsql_my.sh -h juniper.health.unm.edu -n tcrd -c \
	-q "SELECT * FROM expression WHERE etype = 'GTEx'" \
	|gzip -c \
	>$DATADIR/tcrd_expression-gtex.tsv.gz
#
