#!/usr/bin/env Rscript
#
library(readr)
library(data.table)
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==4) {
  (ensgfile <- args[1])
  (hugofile <- args[2])
  (biomartfile <- args[3])
  (ofile <- args[4])
} else if (length(args)==0) {
  ensgfile <- "data/gtex_rnaseq.ensg"
  hugofile <- "data/hugo_protein-coding_gene.tsv"
  biomartfile <- "data/biomart_ENSG2xrefs_human.tsv"
  ofile <- "data/gtex_gene_xref.tsv"
} else {
  message("ERROR: Syntax: gtex_gene_xref.R ENSGFILE HUGOFILE BIOMARTFILE OFILE\n\t...or no args for defaults.")
  quit()
}
message(sprintf("Input ENSG file: %s", ensgfile))
message(sprintf("Input HUGO file: %s", hugofile))
message(sprintf("Input BioMart file: %s", biomartfile))
message(sprintf("Output: %s", ofile))

###
# Ensembl.org/biomart reports 64561 total genes in the Homo sapiens (GRCh38.p12)
# dataset. Selecting only those with Ensembl Protein Family IDs, the number is 22710,
# the protein-encoding genes.
###
biomart <- read_delim(biomartfile, "\t", na=c("", "NA"))
setDT(biomart)
setnames(biomart, c('ENSG', 'ENSGV', 'NCBI', 'HGNCID', 'HGNCSYMB'))
###
# All ENSG (un-versioned) IDs from GTEx RNAseq (1 column).
###
gtex_ensg <- read_delim(ensgfile, col_names=F, delim="\t")
gtex_ensg <- data.table(ENSG = unique(gtex_ensg[[1]]))
n_gtex_ensg <- nrow(gtex_ensg)
message(sprintf("GTEx ENSGs: %d", n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, NCBI)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("ENSGs mapped to NCBI: %d / %d (%.1f%%)",
	sum(!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$NCBI))/n_gtex_ensg))
#
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, HGNCID)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("ENSGs mapped to HGNCID: %d / %d (%.1f%%)",
	sum(!is.na(x$HGNCID)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID))/n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[, .(ENSG, NCBI, HGNCID)]), by='ENSG', all.x=T, all.y=F)
message(sprintf("ENSGs mapped to NCBI and HGNCID: %d / %d (%.1f%%)", 
	sum(!is.na(x$HGNCID)&!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID)&!is.na(x$NCBI))/n_gtex_ensg))
###
# https://www.genenames.org/cgi-bin/statistics
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt
###
hugo <- read_delim(hugofile, "\t")
setDT(hugo)
hugo <- unique(hugo[, .(hgnc_id, symbol, name, gene_family, ensembl_gene_id, uniprot_ids, location)])
setnames(hugo, old=c("uniprot_ids", "location"), new=c("uniprot", "chr"))
#
gtex_gene <- x
###
x <- merge(gtex_gene[!is.na(HGNCID)], hugo, by.x='HGNCID', by.y='hgnc_id', all.x=T, all.y=F)
message(sprintf("ENSGs mapped to HGNC symbols: %d / %d (%.1f%%)",
	sum(!is.na(x$symbol)), n_gtex_ensg,
	100*sum(!is.na(x$symbol))/n_gtex_ensg))
#
message(sprintf("ENSGs mapped to UniProt IDs: %d / %d (%.1f%%)",
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
gtex_gene <- gtex_gene[!(is.na(HGNCID) | is.na(NCBI) | is.na(uniprot))]
gtex_gene <- gtex_gene[ , .(ENSG, NCBI, HGNCID, chr, uniprot, symbol, name)]
message(sprintf("N_ENSG: %d", uniqueN(gtex_gene$ENSG)))
message(sprintf("N_NCBI: %d", uniqueN(gtex_gene$NCBI)))
message(sprintf("N_HGNCID: %d", uniqueN(gtex_gene$HGNCID)))
message(sprintf("N_chromosomal_location: %d", uniqueN(gtex_gene$chr)))
message(sprintf("N_UniProt: %d", uniqueN(gtex_gene$uniprot)))
message(sprintf("N_gene_symbol: %d", uniqueN(gtex_gene$symbol)))
message(sprintf("N_gene_name: %d", uniqueN(gtex_gene$name)))
#
write_delim(gtex_gene, ofile, "\t")
###
