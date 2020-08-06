#!/usr/bin/env Rscript
#
library(readr)
library(data.table)
###
# Ensembl.org/biomart reports 64561 total genes in the Homo sapiens (GRCh38.p12)
# dataset. Selecting only those with Ensembl Protein Family IDs, the number is 22710,
# the protein-encoding genes.
###
biomart <- read_delim("data/ensembl_biomart.tsv", "\t", na=c("","NA"))
setDT(biomart)
setnames(biomart, c('ENSG','ENSGV','HGNCID','NCBI'))
###
# All ENSG (un-versioned) IDs from GTEx RNAseq (1 column).
###
gtex_ensg <- read_delim("data/gtex_rnaseq.ensg", col_names=F, delim="\t")
gtex_ensg <- data.table(ENSG = unique(gtex_ensg[[1]]))
n_gtex_ensg <- nrow(gtex_ensg)
message(sprintf("GTEx ENSGs: %d", n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, NCBI)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("GTEx ENSGs mapped to NCBI: %d / %d (%.1f%%)",
	sum(!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$NCBI))/n_gtex_ensg))
#
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, HGNCID)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("GTEx ENSGs mapped to HGNCID: %d / %d (%.1f%%)",
	sum(!is.na(x$HGNCID)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID))/n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, NCBI, HGNCID)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("GTEx ENSGs mapped to NCBI and HGNCID: %d / %d (%.1f%%)", 
	sum(!is.na(x$HGNCID)&!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID)&!is.na(x$NCBI))/n_gtex_ensg))
###
# https://www.genenames.org/cgi-bin/statistics
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt
###
hugo <- read_delim("data/hugo_protein-coding_gene.tsv", "\t")
setDT(hugo)
hugo <- unique(hugo[, .(hgnc_id, symbol, name, gene_family, ensembl_gene_id, uniprot_ids, location)])
setnames(hugo, old=c("uniprot_ids", "location"), new=c("uniprot", "chr"))
#
gtex_gene <- x
###
x <- merge(gtex_gene[!is.na(HGNCID)], hugo, by.x='HGNCID', by.y='hgnc_id', all.x=T, all.y=F)
message(sprintf("GTEx ENSGs mapped to HGNC symbols: %d / %d (%.1f%%)",
	sum(!is.na(x$symbol)), n_gtex_ensg,
	100*sum(!is.na(x$symbol))/n_gtex_ensg))
#
message(sprintf("GTEx ENSGs mapped to UniProt IDs: %d / %d (%.1f%%)",
                   sum(!is.na(x$uniprot)), n_gtex_ensg,
                   100*sum(!is.na(x$uniprot))/n_gtex_ensg))
#
# ENSG check:
ensg2ensg <- x[!is.na(ensembl_gene_id), .(ENSG, ensembl_gene_id)]
ensg2ensg_mismatch <- ensg2ensg[ENSG!=ensembl_gene_id]
writeLines(sprintf("Ensembl/HGNC ENSG mismatch: %s != %s", ensg2ensg_mismatch$ENSG, ensg2ensg_mismatch$ensembl_gene_id))
###
#
gtex_gene <- x[, ensembl_gene_id := NULL]
#
gtex_gene <- gtex_gene[!(is.na(HGNCID)|is.na(NCBI)|is.na(uniprot))]
gtex_gene <- gtex_gene[ , .(ENSG, NCBI, HGNCID, chr, uniprot, symbol, name)]
ofile <- "data/gtex_gene_xref.tsv"
message(sprintf("Output file, GTEx/Ensembl/HGNC: %s ; nrow: %d", ofile, nrow(gtex_gene)))
message(sprintf("Output file, GTEx/Ensembl/HGNC, ENSG: %d", length(unique(gtex_gene$ENSG))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, NCBI: %d", length(unique(gtex_gene$NCBI))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, HGNCID: %d", length(unique(gtex_gene$HGNCID))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, chromosomal location: %d", length(unique(gtex_gene$chr))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, UniProt: %d", length(unique(gtex_gene$uniprot))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, symbol: %d", length(unique(gtex_gene$symbol))))
message(sprintf("Output file, GTEx/Ensembl/HGNC, name: %d", length(unique(gtex_gene$name))))
#
write_tsv(gtex_gene, ofile)
###
