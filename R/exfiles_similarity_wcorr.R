#!/usr/bin/env Rscript
#############################################################################
### Rationale of weighted Pearson is to unweight smaller, noise-dominated
### expression values.
#############################################################################
### Process input expression profiles (exfiles) and calculate
### weighted correlation for all pairwise combos.
### Input expression profiles expected format 2 ID cols (ENSG,SEX) followed
### by TPM values, one column per tissue.
#############################################################################
### Comparison groups:
###   - F: F vs F
###   - M: M vs M
###   - C: Combined, non-sexed comparisons
### For Combined, profiles are mean of F, M: C = (F+M)/2
#############################################################################
# Conserve memory by writing results directly to file.
#############################################################################
library(readr)
library(wCorr)
#
NiceTime <- function(sec) { sprintf("%d:%02d:%02d",as.integer((as.integer(sec)%%3600)/60),as.integer(sec/3600),as.integer(sec)%%60) }
#
Sys.time()
t0 <- proc.time()
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) { IFILE <- args[1] } else { 
  IFILE <- "data/exfiles_eps.tsv"
}
if (length(args)>1) { OFILE <- args[2] } else { 
  OFILE <- "data/exfiles_eps_WPearson.tsv"
}
if (length(args)>2) { MIN_COR <- as.numeric(args[3]) } else { 
  MIN_COR <- 0.5
}
if (length(args)>3) { MAX_ANTICOR <- as.numeric(args[4]) } else { 
  MAX_ANTICOR <- -0.5
}
#
N_IDCOLS <- 2
#
writeLines(sprintf("INPUT: %s",IFILE))
writeLines(sprintf("OUTPUT: %s",OFILE))
writeLines(sprintf("MIN_COR: %f",MIN_COR))
writeLines(sprintf("MAX_ANTICOR: %f",MAX_ANTICOR))
writeLines(sprintf("N_IDCOLS: %d",N_IDCOLS))
#
fout <- file(OFILE, "w")
writeLines(paste0(c('ENSGA','SEXA','ENSGB','SEXB','wRho'),collapse='\t'), fout)
###
###
wPearson <- function(A,B) { #Vector version (slow)
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
wPearson_mx <- function(A,B) { #Matrix version
  mapply(weightedCorr, as.list(as.data.frame(t(A))), as.list(as.data.frame(t(B))), 
         weights=as.list(as.data.frame(t((A+B)/2))),
         MoreArgs=list(method="Pearson"))
}
###
#
eps <- read_delim(IFILE, "\t", col_types=cols(SEX=col_character()))
#
writeLines(sprintf("Tissue columns: %d",ncol(eps)-N_IDCOLS))
#
eps <- eps[!duplicated(eps[,c("ENSG","SEX")]),]
eps <- eps[order(eps$ENSG),]
#
eps_f <- eps[eps$SEX=="F",]
eps_m <- eps[eps$SEX=="M",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL
###
#F and M gene sets should be same. Not necessary?
#if (!setequal(eps_f$ENSG, eps_m$ENSG)) {
#  eps_f <- eps_f[eps_f$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
#  eps_m <- eps_m[eps_m$ENSG %in% intersect(eps_m$ENSG,eps_m$ENSG),]
#  eps <- eps[eps$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
#}
###
# Compute combined profiles:
eps_c <- aggregate(eps[,!(names(eps) %in% c("ENSG","SEX"))], by=list(ENSG=eps$ENSG), FUN=mean, na.rm=F)
#
# Must have profiles for each gene.
writeLines(sprintf("Expression profiles (F): %d", length(unique(eps_f$ENSG))))
writeLines(sprintf("Expression profiles (M): %d", length(unique(eps_m$ENSG))))
writeLines(sprintf("Expression profiles (C): %d", length(unique(eps_c$ENSG))))
###
#
###
N=nrow(eps_f)
writeLines(sprintf("N = %d ; N_calc_max per group = %d (N(N-1)/2)", N, N*(N-1)/2))
###
n_calc_max <- N*(N-1)/2
#
###
n_calc_total <- 0
n_na_total <- 0
n_ok_total <- 0
###
#F:
eps_mx_f <- as.matrix(eps_f[,N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
#
for (ensgA in eps_f$ENSG) {
  ensgs_this <- eps_f$ENSG[eps_f$ENSG>ensgA]
  group <- "F"
  epAs <- eps_mx_f[rep(ensgA,length(ensgs_this)),]
  epBs <- eps_mx_f[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,group,ensgs_this[results_ok],group,results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", group, n_calc, n_calc_max,100*n_calc/n_calc_max, NiceTime((proc.time()-t0)[3])))
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
###
#M:
eps_mx_m <- as.matrix(eps_m[,N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
#
for (ensgA in eps_m$ENSG) {
  ensgs_this <- eps_m$ENSG[eps_m$ENSG>ensgA]
  group <- "M"
  epAs <- eps_mx_m[rep(ensgA,length(ensgs_this)),]
  epBs <- eps_mx_m[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,group,ensgs_this[results_ok],group,results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", group, n_calc, n_calc_max,100*n_calc/n_calc_max, NiceTime((proc.time()-t0)[3])))
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
#
###
#C:
eps_mx_c <- as.matrix(eps_c[,N_IDCOLS:ncol(eps_c)])
rownames(eps_mx_c) <- eps_c$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
#
for (ensgA in eps_c$ENSG) {
  ensgs_this <- eps_c$ENSG[eps_c$ENSG>ensgA]
  group <- "C"
  epAs <- eps_mx_c[rep(ensgA,length(ensgs_this)),]
  epBs <- eps_mx_c[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,group,ensgs_this[results_ok],group,results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  writeLines(sprintf("Progress (%s): %7d / %7d (%.1f%%) ; elapsed: %s", group, n_calc, n_calc_max,100*n_calc/n_calc_max, NiceTime((proc.time()-t0)[3])))
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
#
###
close(fout)
#
writeLines(sprintf("TOTAL results: calculated: %d ; NAs: %d ; after filtering: %d", n_calc_total, n_na_total, n_ok_total))
#
Sys.time()
#
