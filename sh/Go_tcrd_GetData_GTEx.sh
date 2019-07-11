#!/bin/bash
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
# (Not all columns relevant for GTEx.)
${cwd}/sh/runsql_my.sh -h juniper.health.unm.edu -n tcrd -c \
	-q "SELECT etype, protein_id, tissue, qual_value, number_value, age, sex 
FROM expression WHERE etype = 'GTEx'" \
	|perl -pe 's/\tNULL\t/\t\t/g' \
	|gzip -c \
	>$DATADIR/tcrd_expression-gtex.tsv.gz
#
