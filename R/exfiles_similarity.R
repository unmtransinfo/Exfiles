#!/usr/bin/env Rscript
#############################################################################
### Process input expression profiles (exfiles) and calculate
### one or more measures.
#############################################################################
library(readr)
library(wCorr)
library(dplyr)

###
# Unweight smaller, noise-dominated expression values.
###
wPearson <- function(A,B) { #Vector version (slow)
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
#
wPearson_mx <- function(A,B) { #Matrix version
  mapply(weightedCorr, as.list(as.data.frame(t(A))), as.list(as.data.frame(t(B))), 
         weights=as.list(as.data.frame(t((A+B)/2))),
         MoreArgs=list(method="Pearson"))
}
#
###
eps <- read_delim("data/gtex_rnaseq_prep_profiles.tsv", "\t")
#
genes <- read_delim("data/biomart_ENSG2NCBI.tsv", "\t")
genes <- genes[,c(2,4)]
colnames(genes) <- c('ENSG','HGNC')
eps <- merge(eps, genes, by='ENSG', all.x=T, all.y=F)
eps$ENSG <- eps$HGNC
eps$HGNC <- NULL
colnames(eps)[colnames(eps)=='ENSG'] <- 'gene'

eps <- eps[!duplicated(eps[,c("gene","SEX")]),]
#
#F and M gene sets should be same.
eps_f <- eps[eps$SEX=="F",]
eps_m <- eps[eps$SEX=="M",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL

# Unfortunately must index by gene symbols for now.
if (!setequal(eps_f$gene, eps_m$gene)) {
  eps_f <- eps_f[eps_f$gene %in% intersect(eps_f$gene,eps_m$gene),]
  eps_m <- eps_m[eps_m$gene %in% intersect(eps_m$gene,eps_m$gene),]
  eps <- eps[eps$gene %in% intersect(eps_f$gene,eps_m$gene),]
}
#eps_f <- eps_f[!duplicated(eps_f$gene),]
#eps_m <- eps_m[!duplicated(eps_m$gene),]

###
eps_mx_f <- as.matrix(eps_f[,2:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$gene
eps_mx_m <- as.matrix(eps_m[,2:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$gene
# Populate FM matrix with each named row an expression profile.  F+M/2 mean.
eps_mx_fm <- (as.matrix(eps_f[,2:ncol(eps_f)]) + as.matrix(eps_m[,2:ncol(eps_m)]))/2
rownames(eps_mx_fm) <- eps_f$gene
#
# Must have profiles for each gene in ggc.
writeLines(sprintf("DEBUG: Expression profiles (F): %d", length(unique(eps_f$gene))))
writeLines(sprintf("DEBUG: Expression profiles (M): %d", length(unique(eps_m$gene))))
###
# Gene-gene correlation
gnames <- unique(eps$gene)
gnames <- gnames[order(gnames)]
###
# Conserve memory by writing results directly to file.
###
source(paste0(Sys.getenv("HOME"),"/lib/R/time_utils.R"))
fout <- file("data/exfiles_ggc.csv", "w")
writeLines(paste0(c('Ga','Gb','Cluster','wRho'),collapse=','), fout)
###
t0 <- proc.time()
### F:
i <- 1
n_calc <- 0
n_calc_total <- length(gnames)*(length(gnames)-1)/2
n_na <- 0
cluster <- "F"
for (gA in gnames) {
  i <- i + 1
  gnames_this <- gnames[gnames>gA]
  epAs <- eps_mx_f[rep(gA,length(gnames_this)),]
  epBs <- eps_mx_f[gnames_this,]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s,%s,%s,%.3f",gA,gnames_this,cluster,results_this), fout)
  flush(fout)
  n_calc <- n_calc + length(gnames_this)
  n_na <- n_na + sum(is.na(results_this))
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", cluster, n_calc, n_calc_total,100*n_calc/n_calc_total, time_utils$NiceTime((proc.time()-t0)[3])))
}
writeLines(sprintf("Values calculated: %d ; NAs: %d", n_calc, n_na))
###
### M:
i <- 1
n_calc <- 0
n_calc_total <- length(gnames)*(length(gnames)-1)/2
n_na <- 0
cluster <- "M"
for (gA in gnames) {
  i <- i + 1
  gnames_this <- gnames[gnames>gA]
  epAs <- eps_mx_m[rep(gA,length(gnames_this)),]
  epBs <- eps_mx_m[gnames_this,]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s,%s,%s,%.3f",gA,gnames_this,cluster,results_this), fout)
  flush(fout)
  n_calc <- n_calc + length(gnames_this)
  n_na <- n_na + sum(is.na(results_this))
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", cluster, n_calc, n_calc_total,100*n_calc/n_calc_total, time_utils$NiceTime((proc.time()-t0)[3])))
}
writeLines(sprintf("Values calculated: %d ; NAs: %d", n_calc, n_na))
###
### FM:
i <- 1
n_calc <- 0
n_calc_total <- length(gnames)*(length(gnames)-1)/2
n_na <- 0
cluster <- "FM"
for (gA in gnames) {
  i <- i + 1
  gnames_this <- gnames[gnames>gA]
  epAs <- eps_mx_fm[rep(gA,length(gnames_this)),]
  epBs <- eps_mx_fm[gnames_this,]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s,%s,%s,%.3f",gA,gnames_this,cluster,results_this), fout)
  flush(fout)
  n_calc <- n_calc + length(gnames_this)
  n_na <- n_na + sum(is.na(results_this))
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", cluster, n_calc, n_calc_total,100*n_calc/n_calc_total, time_utils$NiceTime((proc.time()-t0)[3])))
}
writeLines(sprintf("Values calculated: %d ; NAs: %d", n_calc, n_na))
###
close(fout)
###
###
