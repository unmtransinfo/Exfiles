#!/usr/bin/env Rscript
#############################################################################
### Process input expression profiles (exfiles) and calculate
### weighted correlation for all pairwise combos.
#############################################################################
library(readr)
library(wCorr)
#library(dplyr)

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
eps <- eps[!duplicated(eps[,c("ENSG","SEX")]),]
#
#F and M gene sets should be same.
eps_f <- eps[eps$SEX=="female",]
eps_m <- eps[eps$SEX=="male",]
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
eps_mx_f <- as.matrix(eps_f[,2:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
eps_mx_m <- as.matrix(eps_m[,2:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
###
# Conserve memory by writing results directly to file.
###
source(paste0(Sys.getenv("HOME"),"/lib/R/time_utils.R"))
fout <- file("data/gtex_rnaseq_profiles_WPearson.tsv", "w")
writeLines(paste0(c('ENSGA','SEXA','ENSGB','SEXB','wRho'),collapse='\t'), fout)
#
ensgs <- unique(eps$ENSG)
ensgs <- ensgs[order(ensgs)]
###
t0 <- proc.time()
###
i <- 1
n_calc <- 0
n_calc_total <- (2*length(ensgs))*(2*length(ensgs)-1)/2
n_na <- 0
for (ensgA in ensgs) {
  i <- i + 1
  ensgs_this <- ensgs[ensgs>ensgA]
  #FF:
  sexes <- c("F","F")
  epAs <- eps_mx_f[rep(ensgA,length(ensgs_this)),]
  epBs <- eps_mx_f[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this,sexes[2],results_this), fout)
  n_calc <- n_calc + length(ensgs_this)
  n_na <- n_na + sum(is.na(results_this))
  #FM:
  sexes <- c("F","M")
  epBs <- eps_mx_m[ensgs_this,]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this,sexes[2],results_this), fout)
  n_calc <- n_calc + length(ensgs_this)
  n_na <- n_na + sum(is.na(results_this))
  #MM:
  sexes <- c("M","M")
  epAs <- eps_mx_m[rep(ensgA,length(ensgs_this)),]
  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgs_this,sexes[2],results_this), fout)
  n_calc <- n_calc + length(ensgs_this)
  n_na <- n_na + sum(is.na(results_this))
  #
  flush(fout)
  writeLines(sprintf("Progress: %7d / %7d (%.1f%%) ; elapsed: %s", n_calc, n_calc_total,100*n_calc/n_calc_total, time_utils$NiceTime((proc.time()-t0)[3])))
}
writeLines(sprintf("Values calculated: %d ; NAs: %d", n_calc, n_na))
###
close(fout)
###
###
