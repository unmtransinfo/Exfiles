#!/bin/bash
#
cwd=$(pwd)
#
DATADIR="${cwd}/data"
#
(mysql -B -h tcrd.kmc.io -u tcrd -p tcrd660 <<__EOF__
SELECT protein_id, tissue, uberon_id, tpm_level, tpm_level_bysex, tpm_f, tpm_m, tau, tau_bysex, gender AS "sex"
FROM gtex
__EOF__
) \
	|gzip -c \
	>$DATADIR/tcrd_expression-gtex.tsv.gz
#
