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
###   - N: Non-sexed comparisons
### For Non-sexed, profiles N = (F+M)/2
#############################################################################
# Conserve memory by writing results directly to file.
#############################################################################
library(readr)
library(data.table)
library(wCorr)
#
t_start <- Sys.time()
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) { IFILE <- args[1] } else { 
  IFILE <- "data/exfiles_eps.tsv.gz"
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
message(sprintf("INPUT: %s", IFILE))
message(sprintf("OUTPUT: %s", OFILE))
message(sprintf("MIN_COR: %f", MIN_COR))
message(sprintf("MAX_ANTICOR: %f", MAX_ANTICOR))
message(sprintf("N_IDCOLS: %d", N_IDCOLS))
#
fout <- file(OFILE, "w")
writeLines(paste0(c('ENSGA', 'SEXA', 'ENSGB', 'SEXB', 'wRho'), collapse='\t'), fout)
###
###
wPearson <- function(A, B) { #Vector version (slow)
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
wPearson_mx <- function(A, B) { #Matrix version
  mapply(weightedCorr, as.list(as.data.frame(t(A))), as.list(as.data.frame(t(B))), 
         weights=as.list(as.data.frame(t((A+B)/2))),
         MoreArgs=list(method="Pearson"))
}
###
#
eps <- read_delim(IFILE, "\t", col_types=cols(SEX=col_character()))
setDT(eps)
#
message(sprintf("Tissue columns: %d", ncol(eps)-N_IDCOLS))
#
eps <- eps[!duplicated(eps[, .(ENSG, SEX)])]
setorder(eps, ENSG)
#
eps_f <- eps[SEX=="F"]
eps_m <- eps[SEX=="M"]
eps_f[, SEX := NULL]
eps_m[, SEX := NULL]
###
# Compute combined profiles:
eps_c <- eps[, lapply(.SD, mean, na.rm=F), by=ENSG, .SDcols = -c("SEX")]
#
# Must have profiles for each gene.
message(sprintf("Expression profiles (F): %d", uniqueN(eps_f$ENSG)))
message(sprintf("Expression profiles (M): %d", uniqueN(eps_m$ENSG)))
message(sprintf("Expression profiles (N): %d", uniqueN(eps_c$ENSG)))
###
#
###
N <- nrow(eps_f)
message(sprintf("N = %d ; N_calc_max per group = %d (N(N-1)/2)", N, N*(N-1)/2))
###
n_calc_max <- N*(N-1)/2
#
###
n_calc_total <- 0
n_na_total <- 0
n_ok_total <- 0
###
#F:
eps_mx_f <- as.matrix(eps_f[, N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
i <- 0
for (ensgA in eps_f$ENSG) {
  i <- i + 1
  ensgs_this <- eps_f[ENSG>ensgA, ENSG]
  group <- "F"
  epAs <- eps_mx_f[rep(ensgA, length(ensgs_this)), ]
  epBs <- eps_mx_f[ensgs_this, ]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f", ensgA, group, ensgs_this[results_ok], group, results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  if (i %% 100 == 0) {
    t_elapsed <- (Sys.time()-t_start)
    message(sprintf("%s: %7d / %7d (%.1f%%); elapsed: %.1f %s", group, n_calc, n_calc_max, 100*n_calc/n_calc_max, t_elapsed, attr(t_elapsed, "units")))
  }
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
###
#M:
eps_mx_m <- as.matrix(eps_m[, N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
i <- 0
#
for (ensgA in eps_m$ENSG) {
  i <- i + 1
  ensgs_this <- eps_m[ENSG>ensgA, ENSG]
  group <- "M"
  epAs <- eps_mx_m[rep(ensgA, length(ensgs_this)), ]
  epBs <- eps_mx_m[ensgs_this, ]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f", ensgA, group, ensgs_this[results_ok], group, results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  if (i %% 100 == 0) {
    t_elapsed <- (Sys.time()-t_start)
    message(sprintf("%s: %7d / %7d (%.1f%%) ; elapsed: %.1f %s", group, n_calc, n_calc_max,100*n_calc/n_calc_max, t_elapsed, attr(t_elapsed, "units")))
  }
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
#
###
#N:
eps_mx_n <- as.matrix(eps_c[, N_IDCOLS:ncol(eps_c)])
rownames(eps_mx_n) <- eps_c$ENSG
#
n_calc <- 0
n_na <- 0
n_ok <- 0
i <- 0
#
for (ensgA in eps_c$ENSG) {
  i <- i + 1
  ensgs_this <- eps_c[ENSG>ensgA, ENSG]
  group <- "C"
  epAs <- eps_mx_n[rep(ensgA, length(ensgs_this)), ]
  epBs <- eps_mx_n[ensgs_this, ]
  results_this <- wPearson_mx(epAs, epBs)
  n_na <- n_na + sum(is.na(results_this))
  n_calc <- n_calc + length(results_this)
  results_ok <- (!is.na(results_this))&((results_this>=MIN_COR)|(results_this<=MAX_ANTICOR))
  if (sum(results_ok)>0) {
    results_this <- results_this[results_ok]
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f", ensgA, group, ensgs_this[results_ok], group, results_this), fout)
    n_ok <- n_ok + sum(results_ok)
  }
  flush(fout)
  if (i %% 100 == 0) {
    t_elapsed <- (Sys.time()-t_start)
    message(sprintf("%s: %7d / %7d (%.1f%%) ; elapsed: %.1f %s", group, n_calc, n_calc_max, 100*n_calc/n_calc_max, t_elapsed, attr(t_elapsed, "units")))
  }
}
n_calc_total <- n_calc_total + n_calc
n_na_total <- n_na_total + n_na
n_ok_total <- n_ok_total + n_ok
#
###
close(fout)
message(sprintf("OUTPUT file written: %s", OFILE))
#
message(sprintf("TOTAL results: calculated: %d ; NAs: %d ; after filtering: %d", n_calc_total, n_na_total, n_ok_total))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2f %s", t_elapsed, attr(t_elapsed, "units")))
#
