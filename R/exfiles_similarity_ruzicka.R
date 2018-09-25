#!/usr/bin/env Rscript
#############################################################################
### Process input expression profiles (exfiles) and calculate
### similarity for all pairwise combos.
### Input expression profiles expected format 2 ID cols (ENSG,SEX) followed
### by TPM values, one column per tissue.
#############################################################################
library(readr)
library(labdsv)
library(reshape2)
#
Sys.time()
t0 <- proc.time()
#
NiceTime <- function(sec) { sprintf("%d:%02d:%02d",as.integer((as.integer(sec)%%3600)/60),as.integer(sec/3600),as.integer(sec)%%60) }
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) { IFILE <- args[1] } else { 
  IFILE <- "data/exfiles_eps.tsv"
}
if (length(args)>1) { OFILE <- args[2] } else { 
  OFILE <- "data/exfiles_eps_Ruzicka.tsv"
}
if (length(args)>2) { MIN_RUZ <- as.numeric(args[3]) } else { 
  MIN_RUZ <- 0.5
}
N_IDCOLS <- 2
#
writeLines(sprintf("INPUT: %s",IFILE))
writeLines(sprintf("OUTPUT: %s",OFILE))
writeLines(sprintf("MIN_RUZ: %f", MIN_RUZ))
writeLines(sprintf("N_IDCOLS: %d", N_IDCOLS))
#
###
#
eps <- read_delim(IFILE, "\t", col_types=cols(SEX=col_character()))
#
writeLines(sprintf("Tissue columns: %d", ncol(eps)-N_IDCOLS))
#
eps <- eps[!duplicated(eps[,c("ENSG","SEX")]),]
eps <- eps[complete.cases(eps),] #No NAs
#
#F and M gene sets should be same.
eps_f <- eps[eps$SEX=="F",]
eps_m <- eps[eps$SEX=="M",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL

# Index by ENSG
if (!setequal(eps_f$ENSG, eps_m$ENSG)) {
  eps_f <- eps_f[eps_f$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
  eps_m <- eps_m[eps_m$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
  eps <- eps[eps$ENSG %in% intersect(eps_f$ENSG,eps_m$ENSG),]
}
eps <- eps[order(eps$ENSG),]
eps_f <- eps_f[order(eps_f$ENSG),]
eps_m <- eps_m[order(eps_m$ENSG),]
#
# Must have profiles for each gene.
writeLines(sprintf("Expression profiles (F): %d", length(unique(eps_f$ENSG))))
writeLines(sprintf("Expression profiles (M): %d", length(unique(eps_m$ENSG))))
###
#
# To compute FF, MM and FM via one combined matrix.
eps_mx_f <- as.matrix(eps_f[,N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- paste0(eps_f$ENSG, "_F")
eps_mx_m <- as.matrix(eps_m[,N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- paste0(eps_m$ENSG, "_M")
eps_mx_fm <- rbind(eps_mx_f, eps_mx_m)
###
N <- nrow(eps_mx_fm)
writeLines(sprintf("N = %d ; theoretical N_results = %d (N(N-1)/2)", N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_fm, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
#
names(ruz) <- c("ENSG_SEXA", "ENSG_SEXB", "Ruzicka")
ruz$ENSG_SEXA <- as.character(ruz$ENSG_SEXA)
ruz$ENSG_SEXB <- as.character(ruz$ENSG_SEXB)
ruz <- ruz[ruz$ENSG_SEXA<ruz$ENSG_SEXB,] #Effect: select upper matrix.
n_total <- nrow(ruz)
writeLines(sprintf("Results: %d", nrow(ruz)))
ruz <- ruz[ruz$Ruzicka>=MIN_RUZ,]
writeLines(sprintf("Results, post-filtered: %d (%.1f%%)", nrow(ruz),100*nrow(ruz)/n_total))
#
ruz['ENSGA'] <- sub("_.*$", "", ruz$ENSG_SEXA)
ruz['SEXA'] <- sub("^.*_", "", ruz$ENSG_SEXA)
ruz['ENSGB'] <- sub("_.*$", "", ruz$ENSG_SEXB)
ruz['SEXB'] <- sub("^.*_", "", ruz$ENSG_SEXB)
ruz <- ruz[,c("ENSGA","SEXA","ENSGB","SEXB","Ruzicka")]
#
ruz$Ruzicka <- round(ruz$Ruzicka, digits=3)
#
write_delim(ruz, path=OFILE, delim="\t")
###
writeLines(sprintf("Total elapsed: %s", NiceTime((proc.time()-t0)[3])))
Sys.time()
#
