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
eps <- eps[order(eps$ENSG),]
eps <- eps[!duplicated(eps[,c("ENSG","SEX")]),]
eps <- eps[complete.cases(eps),] #No NAs
#
eps_f <- eps[eps$SEX=="F",]
eps_m <- eps[eps$SEX=="M",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL
#
# Compute combined profiles:
eps_c <- aggregate(eps[,!(names(eps) %in% c("ENSG","SEX"))], by=list(ENSG=eps$ENSG), FUN=mean, na.rm=F)
#
writeLines(sprintf("Expression profiles (F): %d", length(unique(eps_f$ENSG))))
writeLines(sprintf("Expression profiles (M): %d", length(unique(eps_m$ENSG))))
writeLines(sprintf("Expression profiles (C): %d", length(unique(eps_c$ENSG))))
###
n_calc_all_total <- 0
n_calc_filtered_total <- 0
#
###
#F:
group <- "F"
eps_mx_f <- as.matrix(eps_f[,N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
#
N <- nrow(eps_f)
writeLines(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_f, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
#
names(ruz) <- c("ENSGA", "ENSGB", "Ruzicka")
ruz$ENSGA <- as.character(ruz$ENSGA)
ruz$ENSGB <- as.character(ruz$ENSGB)
ruz <- ruz[ruz$ENSGA<ruz$ENSGB,] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[ruz$Ruzicka>=MIN_RUZ,]
n_calc_filtered <- nrow(ruz)
writeLines(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz$SEXA <- group
ruz$SEXB <- group
ruz <- ruz[,c("ENSGA","SEXA","ENSGB","SEXB","Ruzicka")]
ruz$Ruzicka <- round(ruz$Ruzicka, digits=3)
write_delim(ruz, path=OFILE, delim="\t")
#
###
#M:
group <- "M"
eps_mx_m <- as.matrix(eps_m[,N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
#
N <- nrow(eps_m)
writeLines(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_m, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
#
names(ruz) <- c("ENSGA", "ENSGB", "Ruzicka")
ruz$ENSGA <- as.character(ruz$ENSGA)
ruz$ENSGB <- as.character(ruz$ENSGB)
ruz <- ruz[ruz$ENSGA<ruz$ENSGB,] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[ruz$Ruzicka>=MIN_RUZ,]
n_calc_filtered <- nrow(ruz)
writeLines(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz$SEXA <- group
ruz$SEXB <- group
ruz <- ruz[,c("ENSGA","SEXA","ENSGB","SEXB","Ruzicka")]
ruz$Ruzicka <- round(ruz$Ruzicka, digits=3)
write_delim(ruz, path=OFILE, delim="\t", append=T)
#
###
#N:
group <- "n"
eps_mx_n <- as.matrix(eps_n[,N_IDCOLS:ncol(eps_n)])
rownames(eps_mx_n) <- eps_n$ENSG
#
###
N <- nrow(eps_n)
writeLines(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_m, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
#
names(ruz) <- c("ENSGA", "ENSGB", "Ruzicka")
ruz$ENSGA <- as.character(ruz$ENSGA)
ruz$ENSGB <- as.character(ruz$ENSGB)
ruz <- ruz[ruz$ENSGA<ruz$ENSGB,] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[ruz$Ruzicka>=MIN_RUZ,]
n_calc_filtered <- nrow(ruz)
writeLines(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz$SEXA <- group
ruz$SEXB <- group
ruz <- ruz[,c("ENSGA","SEXA","ENSGB","SEXB","Ruzicka")]
ruz$Ruzicka <- round(ruz$Ruzicka, digits=3)
write_delim(ruz, path=OFILE, delim="\t", append=T)
###
#
writeLines(sprintf("TOTAL results: all: %d; post-filtered: %d (%.1f%%)", n_calc_all_total, n_calc_filtered_total, 100*n_calc_filtered_total/n_calc_all_total))
writeLines(sprintf("Total elapsed: %s", NiceTime((proc.time()-t0)[3])))
Sys.time()
#
