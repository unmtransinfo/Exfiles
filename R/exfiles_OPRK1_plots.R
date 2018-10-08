#!/usr/bin/env Rscript
###
###
library(readr)
library(dplyr)
library(plotly)

###
#
if (file.exists("R/exfiles/exfiles.Rdata")) {
  message(sprintf("Loading dataset from Rdata..."))
  load("R/exfiles/exfiles.Rdata")
} else {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  ###
  # i, SMTS, SMTSD
  ###
  tissue <- read_delim("R/exfiles/exfiles_tissue_order.tsv", "\t")
  ###
  # ENSG, NCBI, HGNCID, chr, uniprot, symbol, name
  ###
  gene <- read_delim("R/exfiles/gtex_gene_xref.tsv", "\t") # gene attributes
  ###
  idg <- read_delim("R/exfiles/gtex_gene_idg.tsv", "\t") # IDG gene/protein attributes
  idg <- idg[,c("accession", "idgTDL", "idgFamily")]
  colnames(idg) <- c("uniprot", "idgTDL", "idgDTO")
  ###
  # ENSG, SEX, tissue.1, tissue.2, etc.
  ###
  eps <- read_delim("R/exfiles/exfiles_eps.tsv", "\t", col_types=cols(SEX=col_character()))
  ###
  # ENSGA, ENSGB, Group, wRho, Ruzicka
  ###
  ggc <- read_delim("R/exfiles/exfiles_ggc.tsv.gz", "\t", col_types="cccdd")
  ggc$Group[ggc$Group=="C"] <- "N"
  ###
  #save(tissue, gene, idg, eps, ggc, file="R/exfiles/exfiles.Rdata")
}
#
tissue_missing <- setdiff(tissue$SMTSD, colnames(eps))
if (length(tissue_missing)>0) {
  message(sprintf("NOTE: TISSUE_MISSING: %d. %s\n", 1:length(tissue_missing),
tissue_missing
))
} else {
  message(sprintf("All tissues found."))
}
tissue <- tissue[tissue$SMTSD %in% colnames(eps),]
eps <- eps[,c("ENSG","SEX",tissue$SMTSD)]
#
message(sprintf("Tissue count: %d",nrow(tissue)))
#
ensgs <- intersect(eps$ENSG, c(ggc$ENSGA,ggc$ENSGB))
ensgs <- intersect(ensgs, gene$ENSG)
message(sprintf("Gene count: %d",length(ensgs)))
#
ensg_dups <- gene$ENSG[duplicated(gene$ENSG)]
message(sprintf("Duplicated/ambiguous gene IDs: %d", length(ensg_dups)))
#
gene <- gene[gene$ENSG %in% ensgs,]
#
message(sprintf("Unknown/unmapped gene SYMBs: %d", sum(is.na(gene$symbol))))
#
gene <- gene[!is.na(gene$symbol),]
#
symb_dups <- gene$symbol[duplicated(gene$symbol)]
message(sprintf("Duplicated/ambiguous gene SYMBs: %d", length(symb_dups)))
#
gene <- gene[!duplicated(gene$ENSG),]
gene <- gene[!duplicated(gene$symbol),]
#
message(sprintf("Gene count: %d", nrow(gene)))
message(sprintf("Gene unique ENSG count: %d", length(unique(gene$ENSG))))
message(sprintf("Gene unique SYMB count: %d", length(unique(gene$symbol))))
#
gene <- merge(gene, idg, by="uniprot", all.x=T, all.y=F)
message(sprintf("Gene unique UniProt count: %d", length(unique(gene$uniprot))))
#
tbl <- table(gene$idgDTO)
message(sprintf("%24s: %6s\n", names(tbl), tbl))
#
#############################################################################
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  wCorr::weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
Ruzicka <- function(A,B) {
  sum(pmin(A,B), rm.na=T)/sum(pmax(A,B), rm.na=T)
}
###
EpLogIf <- function(ep, condition) {
  if (condition) { return(log10(ep+1)) }
  else { return(ep) }
}
#############################################################################
###
qryA <- "OPRK1"
ensgA <- gene$ENSG[gene$symbol==qryA]
#
LogY <- T
#
xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10), categoryorder = "array", categoryarray = tissue$SMTSD)
yaxis = list(title="Expression: LOG<SUB>10</SUB>(1+TPM)")
###
###
qryA_profile_f <- as.numeric(eps[eps$ENSG==ensgA & eps$SEX=="F",][1,tissue$SMTSD])
qryA_profile_m <- as.numeric(eps[eps$ENSG==ensgA & eps$SEX=="M",][1,tissue$SMTSD])
qryA_profile_n <- (qryA_profile_f+qryA_profile_m)/2
###
#
###
# ggc includes only high correlations, so more reliable to compute from profiles.
###
### GeneB query:
qryB_regex <- "^ADCY[0-9]"
#
gene_hits <- gene[grepl(qryB_regex, gene$symbol),]
gene_hits <- gene_hits[order(as.numeric(sub("ADCY", "", gene_hits$symbol))),]
#
###
#N:
p_n <- plot_ly() %>%
  add_trace(name = paste("(N)", qryA), x = tissue$SMTSD, y = EpLogIf(qryA_profile_n, LogY),  
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol = "circle", size = 10),
            text = paste0(qryA, ": ", gene$name[gene$ENSG==ensgA])
  ) %>%
  layout(xaxis = xaxis, yaxis = yaxis, 
        title = sprintf("GTEx Gene-Tissue Profiles (N): %s vs %s", qryA, qryB_regex),
        margin = list(t=100,r=80,b=160,l=60),
        legend = list(x = .9, y = 1),
        font = list(family = "Arial", size=14))
#
p_f <- plot_ly() %>%
  add_trace(name = paste("(F)", qryA), x = tissue$SMTSD, y = EpLogIf(qryA_profile_f, LogY),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol = "circle", size = 10),
            text = paste0(qryA, ": ", gene$name[gene$ENSG==ensgA])
  ) %>%
  layout(xaxis = xaxis, yaxis = yaxis, 
        title = sprintf("GTEx Gene-Tissue Profiles (F): %s vs %s", qryA, qryB_regex),
        margin = list(t=100,r=80,b=160,l=60),
        legend = list(x = .9, y = 1),
        font = list(family = "Arial", size=14))
#
p_m <- plot_ly() %>%
  add_trace(name = paste("(M)", qryA), x = tissue$SMTSD, y = EpLogIf(qryA_profile_m, LogY),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol = "circle", size = 10),
            text = paste0(qryA, ": ", gene$name[gene$ENSG==ensgA])
  ) %>%
  layout(xaxis = xaxis, yaxis = yaxis, 
         title = sprintf("GTEx Gene-Tissue Profiles (M): %s vs %s", qryA, qryB_regex),
         margin = list(t=100,r=80,b=160,l=60),
         legend = list(x = .9, y = 1),
         font = list(family = "Arial", size=14))
#
for (i in 1:nrow(gene_hits)) {
  hit <- gene_hits$ENSG[i]
  hit_symbol <- gene_hits$symbol[i]
  hit_name <- gene_hits$name[i]
  hit_profile_f <- as.numeric(eps[eps$ENSG==hit & eps$SEX=="F",][1,tissue$SMTSD])
  hit_profile_m <- as.numeric(eps[eps$ENSG==hit & eps$SEX=="M",][1,tissue$SMTSD])
  hit_profile_n <- (hit_profile_f+hit_profile_m)/2
  rhoNah <- wPearson(qryA_profile_n, hit_profile_n)
  ruzNah <- Ruzicka(qryA_profile_n, hit_profile_n)
  rhoFah <- wPearson(qryA_profile_f, hit_profile_f)
  ruzFah <- Ruzicka(qryA_profile_f, hit_profile_f)
  rhoMah <- wPearson(qryA_profile_m, hit_profile_m)
  ruzMah <- Ruzicka(qryA_profile_m, hit_profile_m)
  writeLines(sprintf("hit: %s; %s; rhoN=%f; ruzN=%f", hit, hit_symbol, rhoNah, ruzNah))
  writeLines(sprintf("hit: %s; %s; rhoF=%f; ruzF=%f", hit, hit_symbol, rhoFah, ruzFah))
  writeLines(sprintf("hit: %s; %s; rhoM=%f; ruzM=%f", hit, hit_symbol, rhoMah, ruzMah))
  p_n <- add_trace(p_n, name=paste("(N)", hit_symbol), x=tissue$SMTSD, y=EpLogIf(hit_profile_n, LogY),
            type='scatter', mode='lines+markers',
            marker=list(symbol="circle", size=10),
            text=paste0(hit_symbol, ": ", hit_name))
  p_f <- add_trace(p_f, name=paste("(F)", hit_symbol), x=tissue$SMTSD, y=EpLogIf(hit_profile_f, LogY),
            type='scatter', mode='lines+markers',
            marker=list(symbol="circle", size=10),
            text=paste0(hit_symbol, ": ", hit_name))
  p_m <- add_trace(p_m, name=paste("(F)", hit_symbol), x=tissue$SMTSD, y=EpLogIf(hit_profile_m, LogY),
            type='scatter', mode='lines+markers',
            marker=list(symbol="circle", size=10),
            text=paste0(hit_symbol, ": ", hit_name))
}
#
p_n
p_f
p_m
#
###

