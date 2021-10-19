#!/usr/bin/env Rscript
#######################################################################################
### Inputs:
###	exfiles_tissue_order.tsv - manually curated, from GDoc
###	gtex_gene_xref.tsv - from gtex_gene_xref.R
###	tcrd_targets.tsv - BioClients.idg.Client 
###	exfiles_eps.tsv - expression profiles; gtex_rnaseq_prep_app.py
###	exfiles_ggc.tsv - gene-gene comparisons; exfiles_similarity_post.py
### Output:
###	exfiles.Rdata - For Shiny app.py.
#######################################################################################
library(readr)
library(data.table, quietly=T)
#
###
# This code runs once for all sessions.
###
GTEX_RELEASE <- "v8 (2017)"
#
###
t0 <- proc.time()
message(sprintf("Loading dataset from files, writing Rdata..."))
###
# i, SMTS, SMTSD
###
tissue <- read_delim("data/exfiles_tissue_order.tsv", "\t", col_types=cols(.default=col_character(), SEX_SPECIFIC=col_logical()))
setDT(tissue)
###
# ENSG, NCBI, HGNCID, chr, uniprot, symbol, name
###
gene <- read_delim("data/gtex_gene_xref.tsv", "\t") # gene attributes
setDT(gene)
###
idg <- read_delim("data/tcrd_targets.tsv", "\t") # IDG gene attributes
setDT(idg)
idg <- idg[, .(uniprotId, TDL, tcrdTargetFamily)]
idg <- idg[!duplicated(uniprotId)]
setnames(idg, c("uniprot", "idgTDL", "idgFamily"))
###
# ENSG, SEX, tissue.1, tissue.2, etc.
###
eps <- read_delim("data/exfiles_eps.tsv", "\t", col_types=cols(.default=col_double(), ENSG=col_character(), SEX=col_character()))
setDT(eps)
###
# Filter tissues absent from custom tissues list.
hiddenTissues <- setdiff(names(eps)[3:ncol(eps)], tissue$SMTSD)
for (tis in hiddenTissues) {
  message(sprintf("Tissue hidden (SMTSD): %s", tis))
}
tissue <- tissue[SMTSD %in% colnames(eps)]
TAGS_THIS <- c("ENSG", "SEX", tissue$SMTSD)
eps <- eps[, ..TAGS_THIS]
###
# ENSGA, ENSGB, Group, wRho, Ruzicka
###
ggc <- read_delim("data/exfiles_ggc.tsv", "\t", col_types="cccdd")
setDT(ggc)
ggc[, Combo := round(wRho*Ruzicka, digits=2)]
#
ensgs <- intersect(eps$ENSG, c(ggc$ENSGA, ggc$ENSGB))
ensgs <- intersect(ensgs, gene$ENSG)
eps <- eps[ENSG %in% ensgs]
ggc <- ggc[(ENSGA %in% ensgs) & (ENSGB %in% ensgs)]
#
gene <- gene[ENSG %in% ensgs]
gene <- merge(gene, idg, by="uniprot", all.x=F, all.y=F)
gene <- gene[!is.na(symbol)]
gene <- gene[!duplicated(ENSG)]
gene <- gene[!duplicated(symbol)]
#
gene_menu <- gene
gene_menu[, symbol := ifelse(!is.na(symbol), symbol, ENSG)] #NAs break autocomplete. 
###
save(tissue, gene, gene_menu, idg, eps, ggc, file="data/exfiles.Rdata")
#
message(sprintf("Gene count (ENSG): %d", uniqueN(gene$ENSG)))
message(sprintf("Gene count (SYMB): %d", uniqueN(gene$symbol)))
message(sprintf("Gene count (UniProt): %d", uniqueN(gene$uniprot)))
message(sprintf("Tissue count (profiles): %d", ncol(eps)-2))
message(sprintf("Tissue count (shown): %d", uniqueN(tissue$SMTSD)))
#
