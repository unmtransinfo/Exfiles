#!/usr/bin/env Rscript
#
library(readr)
###
# Ensembl.org/biomart reports 64561 total genes in the Homo sapiens (GRCh38.p12)
# dataset. Selecting only those with Ensembl Protein Family IDs, the number is 22710,
# the protein-encoding genes.
###
biomart <- read_delim("data/ensembl_biomart.tsv", "\t", na=c("","NA"))
colnames(biomart) <- c('ENSG','ENSGV','HGNCID','NCBI')
###
# All ENSG (un-versioned) IDs from GTEx RNAseq (1 column).
###
gtex_ensg <- read_csv("data/gtex_rnaseq.ensg",col_names=F)
gtex_ensg <- unique(gtex_ensg)
colnames(gtex_ensg) <- 'ENSG'
n_gtex_ensg <- nrow(gtex_ensg)
writeLines(sprintf("GTEx ENSGs: %d", n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[,c('ENSG','NCBI')]), by='ENSG', all.x=T, all.y=F)
writeLines(sprintf("GTEx ENSGs mapped to NCBI: %d / %d (%.1f%%)",
	sum(!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$NCBI))/n_gtex_ensg))
#
x <- merge(gtex_ensg, unique(biomart[,c('ENSG','HGNCID')]), by='ENSG', all.x=T, all.y=F)
writeLines(sprintf("GTEx ENSGs mapped to HGNCID: %d / %d (%.1f%%)",
	sum(!is.na(x$HGNCID)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID))/n_gtex_ensg))
###
x <- merge(gtex_ensg, unique(biomart[,c('ENSG','NCBI','HGNCID')]), by='ENSG', all.x=T, all.y=F)
writeLines(sprintf("GTEx ENSGs mapped to NCBI and HGNCID: %d / %d (%.1f%%)", 
	sum(!is.na(x$HGNCID)&!is.na(x$NCBI)), n_gtex_ensg,
	100*sum(!is.na(x$HGNCID)&!is.na(x$NCBI))/n_gtex_ensg))
###
# https://www.genenames.org/cgi-bin/statistics
# ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt
###
hugo <- read_delim("data/hugo_protein-coding_gene.txt", "\t")
hugo <- unique(hugo[,c("hgnc_id","symbol","name","gene_family","ensembl_gene_id")])
#
gtex_gene <- x
###
x <- merge(gtex_gene[!is.na(x$HGNCID),], hugo, by.x='HGNCID', by.y='hgnc_id', all.x=T, all.y=F)
writeLines(sprintf("GTEx ENSGs mapped to HGNC symbols: %d / %d (%.1f%%)",
	sum(!is.na(x$symbol)), n_gtex_ensg,
	100*sum(!is.na(x$symbol))/n_gtex_ensg))
#
# ENSG check:
ensg2ensg <- x[!is.na(x$ensembl_gene_id),c("ENSG","ensembl_gene_id")]
x$ensembl_gene_id <- NULL
ensg2ensg_mismatch <- ensg2ensg[ensg2ensg$ENSG!=ensg2ensg$ensembl_gene_id,]
writeLines(sprintf("Ensembl/HGNC ENSG mismatch: %s != %s", ensg2ensg_mismatch$ENSG, ensg2ensg_mismatch$ensembl_gene_id))
###
#
gtex_gene <- x
#
gtex_gene <- gtex_gene[!(is.na(x$HGNCID)&is.na(x$NCBI)),]
gtex_gene <- gtex_gene[,c("ENSG","NCBI","HGNCID","symbol","name")]
ofile <- "data/gtex_gene_xref.tsv"
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC: %s ; nrow: %d", ofile, nrow(gtex_gene)))
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC, ENSG: %d", length(unique(gtex_gene$ENSG))))
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC, NCBI: %d", length(unique(gtex_gene$NCBI))))
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC, HGNCID: %d", length(unique(gtex_gene$HGNCID))))
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC, symbol: %d", length(unique(gtex_gene$symbol))))
writeLines(sprintf("Output file, GTEx/Ensembl/HGNC, name: %d", length(unique(gtex_gene$name))))
write_tsv(gtex_gene, ofile)
###
