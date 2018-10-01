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
#
#F and M gene sets should be same.
eps_f <- eps[eps$SEX=="F",]
eps_m <- eps[eps$SEX=="M",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL

# Index by ENSG
if (!setequal(eps_f$ENSG, eps_m$ENSG)) {
  eps_f <- eps_f[eps_f$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
  eps_m <- eps_m[eps_m$ENSG %in% intersect(eps_m$ENSG,eps_m$ENSG),]
  eps <- eps[eps$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
}
# Must have profiles for each gene.
writeLines(sprintf("Expression profiles (F): %d", length(unique(eps_f$ENSG))))
writeLines(sprintf("Expression profiles (M): %d", length(unique(eps_m$ENSG))))
###

### Matrices efficient, faster.
###
eps_mx_f <- as.matrix(eps_f[,N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
eps_mx_m <- as.matrix(eps_m[,N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
###
# Conserve memory by writing results directly to file.
###
#
ensgs <- unique(eps$ENSG)
ensgs <- ensgs[order(ensgs)]
###
N=2*length(ensgs)
writeLines(sprintf("N = %d ; theoretical N_results = %d (N(N-1)/2)", N, N*(N-1)/2))
###
n_calc <- 0
n_calc_total <- N*(N-1)/2
n_na <- 0
n_ok <- 0
for (ensgA in ensgs) {
  ensgs_this <- ensgs[ensgs>ensgA]
  #FF:
  sexes <- c("F","F")
  epAs <- eps_mx_f[rep(ensgA,length(ensgs_this)),]
  epBs <- eps_mx_f[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this[results_ok],sexes[2],results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  #FM:
  sexes <- c("F","M")
  epBs <- eps_mx_m[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this[results_ok],sexes[2],results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  #MM:
  sexes <- c("M","M")
  epAs <- eps_mx_m[rep(ensgA,length(ensgs_this)),]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this[results_ok],sexes[2],results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  #
  flush(fout)
  writeLines(sprintf("Progress: %7d / %7d (%.1f%%) ; elapsed: %s", n_calc, n_calc_total,100*n_calc/n_calc_total, NiceTime((proc.time()-t0)[3])))
}
#
close(fout)
#
writeLines(sprintf("Results calculated: %d ; NAs: %d ; after filtering: %d", n_calc, n_na, n_ok))
#
Sys.time()
#
