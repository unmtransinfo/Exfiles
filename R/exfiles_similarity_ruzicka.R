#!/usr/bin/env Rscript
#############################################################################
### Process input expression profiles (exfiles) and calculate
### similarity for all pairwise combos.
### Input expression profiles expected format 2 ID cols (ENSG,SEX) followed
### by TPM values, one column per tissue.
#############################################################################
library(readr)
library(data.table)
library(reshape2)
library(labdsv)
#
Sys.time()
t_start <- Sys.time()
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
message(sprintf("INPUT: %s",IFILE))
message(sprintf("OUTPUT: %s",OFILE))
message(sprintf("MIN_RUZ: %f", MIN_RUZ))
message(sprintf("N_IDCOLS: %d", N_IDCOLS))
#
###
#
eps <- read_delim(IFILE, "\t", col_types=cols(.default=col_double(), ENSG=col_character(), SEX=col_character()))
setDT(eps)
#
message(sprintf("Tissue columns: %d", ncol(eps)-N_IDCOLS))
#
setorder(eps, ENSG)
eps <- eps[!duplicated(eps[, .(ENSG, SEX)])]
for (j in 3:ncol(eps))
  set(eps, which(is.na(eps[[j]])), j, 0)
#
eps_f <- eps[SEX=="F"]
eps_m <- eps[SEX=="M"]
eps_f[, SEX := NULL]
eps_m[, SEX := NULL]
#
# Compute combined profiles:
eps_c <- eps[, lapply(.SD, mean, na.rm=F), by=ENSG, .SDcols = -c("SEX")] 
#
message(sprintf("Expression profiles (F): %d", uniqueN(eps_f$ENSG)))
message(sprintf("Expression profiles (M): %d", uniqueN(eps_m$ENSG)))
message(sprintf("Expression profiles (C): %d", uniqueN(eps_c$ENSG)))
###
n_calc_all_total <- 0
n_calc_filtered_total <- 0
#
###
#F:
group <- "F"
eps_mx_f <- as.matrix(eps_f[, N_IDCOLS:ncol(eps_f)])
rownames(eps_mx_f) <- eps_f$ENSG
#
N <- nrow(eps_f)
message(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_f, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
setDT(ruz)
#
setnames(ruz, c("ENSGA", "ENSGB", "Ruzicka"))
ruz[, ENSGA := as.character(ENSGA)]
ruz[, ENSGB := as.character(ENSGB)]
ruz <- ruz[ENSGA<ENSGB] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[Ruzicka>=MIN_RUZ]
n_calc_filtered <- nrow(ruz)
message(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz[, SEXA := group]
ruz[, SEXB := group]
ruz <- ruz[, .(ENSGA, SEXA, ENSGB, SEXB, Ruzicka)]
ruz[, Ruzicka := round(Ruzicka, digits=3)]
write_delim(ruz, path=OFILE, delim="\t")
#
###
#M:
group <- "M"
eps_mx_m <- as.matrix(eps_m[, N_IDCOLS:ncol(eps_m)])
rownames(eps_mx_m) <- eps_m$ENSG
#
N <- nrow(eps_m)
message(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_m, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
setDT(ruz)
#
setnames(ruz, c("ENSGA", "ENSGB", "Ruzicka"))
ruz[, ENSGA := as.character(ENSGA)]
ruz[, ENSGB := as.character(ENSGB)]
ruz <- ruz[ENSGA<ENSGB] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[Ruzicka>=MIN_RUZ]
n_calc_filtered <- nrow(ruz)
message(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz$SEXA <- group
ruz$SEXB <- group
ruz <- ruz[, .(ENSGA, SEXA, ENSGB, SEXB, Ruzicka)]
ruz[, Ruzicka := round(Ruzicka, digits=3)]
write_delim(ruz, path=OFILE, delim="\t", append=T)
#
###
#C:
group <- "C"
eps_mx_c <- as.matrix(eps_c[, N_IDCOLS:ncol(eps_c)])
rownames(eps_mx_c) <- eps_c$ENSG
#
###
N <- nrow(eps_c)
message(sprintf("(%s) N = %d ; per group N_calc_max = %d (N(N-1)/2)", group, N, N*(N-1)/2))
#
# Generates full matrix, but only want upper.
ruzd <- labdsv::dsvdis(eps_mx_c, "ruzicka", upper=T) #dist
ruzm <- -as.matrix(ruzd) + 1 #matrix: sim=(1-dist)
ruz <- reshape2::melt(ruzm)
setDT(ruz)
#
setnames(ruz, c("ENSGA", "ENSGB", "Ruzicka"))
ruz[, ENSGA := as.character(ENSGA)]
ruz[, ENSGB := as.character(ENSGB)]
ruz <- ruz[ENSGA<ENSGB] #Effect: select upper matrix.
n_calc_all <- nrow(ruz)
ruz <- ruz[Ruzicka>=MIN_RUZ]
n_calc_filtered <- nrow(ruz)
message(sprintf("Results (%s): all: %d; post-filtered: %d (%.1f%%)", group, n_calc_all, n_calc_filtered, 100*n_calc_filtered/n_calc_all))
n_calc_all_total <- n_calc_all_total + n_calc_all
n_calc_filtered_total <- n_calc_filtered_total + n_calc_filtered
#
ruz$SEXA <- group
ruz$SEXB <- group
ruz <- ruz[, .(ENSGA, SEXA, ENSGB, SEXB, Ruzicka)]
ruz[, Ruzicka := round(Ruzicka, digits=3)]
write_delim(ruz, path=OFILE, delim="\t", append=T)
###
#
message(sprintf("TOTAL results: all: %d; post-filtered: %d (%.1f%%)", n_calc_all_total, n_calc_filtered_total, 100*n_calc_filtered_total/n_calc_all_total))
#
t_elapsed <- (Sys.time()-t_start)
message(sprintf("Elapsed time: %.2f %s", t_elapsed, attr(t_elapsed, "units")))
Sys.time()
