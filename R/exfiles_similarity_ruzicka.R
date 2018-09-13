#!/usr/bin/env Rscript
#############################################################################
#############################################################################
### Process input expression profiles (exfiles) and calculate
### similarity for all pairwise combos.
### Input expression profiles expected format 2 ID cols (ENSG,SEX) followed
### by TPM values, one column per tissue.
#############################################################################
library(readr)
library(labdsv)
#
source(paste0(Sys.getenv("HOME"),"/lib/R/time_utils.R"))
#
t0 <- proc.time()
#
args <- commandArgs(trailingOnly=TRUE)
if (length(args)>0) { IFILE <- args[1] } else { 
  IFILE <- "data/gtex_rnaseq_prep_profiles.tsv"
}
if (length(args)>1) { OFILE <- args[2] } else { 
  OFILE <- "data/gtex_rnaseq_profiles_Ruzicka.tsv"
}
writeLines(sprintf("INPUT: %s",IFILE))
writeLines(sprintf("OUTPUT: %s",OFILE))
#
###
fout <- file(OFILE, "w")
writeLines(paste0(c('ENSGA','SEXA','ENSGB','SEXB','Ruzicka'),collapse='\t'), fout)
###
#
eps <- read_delim(IFILE, "\t")
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
  eps_m <- eps_m[eps_m$ENSG %in% intersect(eps_m$ENSG,eps_m$ENSG),]
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
eps_mx_f <- as.matrix(eps_f[,2:ncol(eps_f)])
rownames(eps_mx_f) <- paste0(eps_f$ENSG, "_F")
eps_mx_m <- as.matrix(eps_m[,2:ncol(eps_m)])
rownames(eps_mx_m) <- paste0(eps_m$ENSG, "_M")
eps_mx_fm <- rbind(eps_mx_f, eps_mx_m)
###
#n_calc_max <- (2*length(ensgs))*(2*length(ensgs)-1)/2
#
n_out <- 0
n_na <- 0
#
ruzd <- labdsv::dsvdis(eps_mx_fm, "ruzicka", upper=T) #dist
ruzs <- -as.matrix(ruzd) + 1 #sim = (1 - dist)
#
for (ensgA_sexA in rownames(ruzs)) {
  for (ensgB_sexB in colnames(ruzs)) {
    if (ensgA_sexA>=ensgB_sexB) { next }
    val <- ruzs[ensgA_sexA,ensgB_sexB]
    if (is.na(val)) {
      n_na <- n_na + 1
      next
    }
    ensgA <- sub("_.*$", "", ensgA_sexA)
    ensgB <- sub("_.*$", "", ensgB_sexB)
    sexA <- sub("^.*_", "", ensgA_sexA)
    sexB <- sub("^.*_", "", ensgB_sexB)
    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexA,ensgB,sexB,val), fout)
    n_out <- n_out + 1
  }
}
#############################################################################
#sexes <- c("F","F")
#ruzd <- labdsv::dsvdis(eps_mx_f, "ruzicka") #dist
#ruzs <- -as.matrix(ruzd) + 1 #sim = (1 - dist)
##
#for (ensgA in rownames(ruzs)) {
#  for (ensgB in colnames(ruzs)) {
#    if (ensgA>=ensgB) { next }
#    val <- ruzs[ensgA,ensgB]
#    if (is.na(val)) {
#      n_na <- n_na + 1
#      next
#    }
#    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgB,sexes[2],val), fout)
#    n_out <- n_out + 1
#  }
#}
##
#sexes <- c("M","M")
#ruzd <- labdsv::dsvdis(eps_mx_m, "ruzicka") #dist
#ruzs <- -as.matrix(ruzd) + 1 #sim = (1 - dist)
##
#for (ensgA in rownames(ruzs)) {
#  for (ensgB in colnames(ruzs)) {
#    if (ensgA>=ensgB) { next }
#    val <- ruzs[ensgA,ensgB]
#    if (is.na(val)) {
#      n_na <- n_na + 1
#      next
#    }
#    writeLines(sprintf("%s\t%s\t%s\t%s\t%.3f",ensgA,sexes[1],ensgB,sexes[2],val), fout)
#    n_out <- n_out + 1
#  }
#}
#############################################################################
###
close(fout)
#
writeLines(sprintf("Values out: %d ; NAs: %d", n_out, n_na))
