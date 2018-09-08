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
         MoreArgs = list( method="Pearson"))
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

#
#F and M gene sets should be same.
eps_f <- eps[eps$SEX=="female",]
eps_m <- eps[eps$SEX=="male",]
eps_f$SEX <- NULL
eps_m$SEX <- NULL

# Unfortunately must index by gene symbols for now.
if (!setequal(eps_f$gene, eps_m$gene)) {
  eps_f <- eps_f[eps_f$gene %in% intersect(eps_f$gene,eps_m$gene),]
  eps_m <- eps_m[eps_m$gene %in% intersect(eps_m$gene,eps_m$gene),]
  eps <- eps[eps$gene %in% intersect(eps_f$gene,eps_m$gene),]
}
eps <- eps[!duplicated(eps[,c("gene","SEX")]),]
eps_f <- eps_f[!duplicated(eps_f$gene),]
eps_m <- eps_m[!duplicated(eps_m$gene),]

###
# Populate matrix with each named row an expression profile.  M_F mean.
eps_mx <- (as.matrix(eps_f[,2:ncol(eps_f)]) + as.matrix(eps_m[,2:ncol(eps_m)]))/2
#M+F mean
rownames(eps_mx) <- eps_f$gene
#
# Must have profiles for each gene in ggc.
writeLines(sprintf("DEBUG: Expression profiles (F): %d",
length(unique(eps_f$gene))))
writeLines(sprintf("DEBUG: Expression profiles (M): %d",
length(unique(eps_m$gene))))
writeLines(sprintf("DEBUG: Expression profiles matrix ((F+M)/2): %d", nrow(eps_mx)))
###
# Gene-gene correlation
gnames <- unique(eps$gene)

ggc <- data.frame(Ga = rep(NA, length(gnames)*length(gnames)/2), Gb = NA, wRho = NA)
i <- 0
for (gA in gnames) {
  for (gB in gnames) {
    if (gA >= gB) { next }
    i <- i + 1
    ggc$Ga[i] <- gA
    ggc$Gb[i] <- gB
  }
}
ggc <- ggc[!is.na(ggc$Ga),]
###
# Conserve memory by writing results directly to file.
###
source(paste0(Sys.getenv("HOME"),"/lib/R/time_utils.R"))
fout <- file("data/exfiles_ggc.csv", "w")
writeLines(paste0(colnames(ggc),collapse=','), fout)
n_chunk <- 1e2
t0 <- proc.time()
i <- 1
n_calc <- 0
n_na <- 0
while (i<nrow(ggc)) {
  i_next <- min(i+n_chunk, nrow(ggc)+1)
  gAs <- ggc$Ga[i:(i_next-1)]
  gBs <- ggc$Gb[i:(i_next-1)]
  epAs <- eps_mx[gAs,]
  epBs <- eps_mx[gBs,]

  results_this <- wPearson_mx(epAs, epBs)
  writeLines(sprintf("%s,%s,%.3f",ggc$Ga[i:(i_next-1)],ggc$Gb[i:(i_next-1)],results_this), fout)
  #ggc$wRho[i:(i_next-1)] <- results_this
  flush(fout)

  n_calc <- n_calc + (i_next-i)
  n_na <- n_na + sum(is.na(results_this))
  
  if ((i %% n_chunk)==1) {
    writeLines(sprintf("Progress: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc),100*i/nrow(ggc), time_utils$NiceTime((proc.time()-t0)[3])))
  }
  i <- i_next
}
writeLines(sprintf("Total: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc), 100*i/nrow(ggc), time_utils$NiceTime((proc.time()-t0)[3])))
###
writeLines(sprintf("Values calculated: %d ; NAs: %d", n_calc, n_na))
###
close(fout)
###
#ggc <- ggc[!is.na(ggc$wRho),]
#ggc$wRho <- round(ggc$wRho, digits=3)
###
#gzout <- gzfile("data/exfiles_ggc.csv.gz","w")
#write.csv(ggc, gzout, row.names=F)
#close(gzout)
###
