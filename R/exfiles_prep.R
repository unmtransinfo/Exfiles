#!/usr/bin/env Rscript
###
### This processes the original files provided by Giovanni (early Feb 2018): 
### (1) result_recall.csv: gene-gene associations
### (2) receptors_all.csv: expression profiles
###
library(readr)
library(wCorr)
library(dplyr)
source(paste0(Sys.getenv("HOME"),"/lib/R/hf_utils.R"))
###
# VECTORIZE FUNCTIONS!
###
RMSD <- function(A,B) {
  sqrt(sum((A-B)^2)/length(A))
}
#
RMSD_mx <- function(A,B) {
  sqrt(rowSums((A-B)^2)/ncol(A))
}
#
###
# Note A%*%B is dot product of two vectors => scalar.
# Note rowSums(A*B) is pairwise dot products of rows.
###
Tanimoto <- function(A,B) { #Vector version (slow)
  (A%*%B) / (A%*%A + B%*%B - A%*%B)
}
#
Tanimoto_mx <- function(A, B) { #Matrix version
  rowSums(A*B) / (rowSums(A*A) + rowSums(B*B) - rowSums(A*B))
}
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
#
### Expresion profiles from Gio.
### LOG10(1+TPM)
eps <- read_delim("data/receptors_all.csv.gz", "\t", col_names=F)
eps <- transform(eps, tissue=do.call(rbind, strsplit(X5, '[, ]+', fixed=F)), stringsAsFactors=F)
eps <- rename(eps, ensg.ver=`X1`, gene=`X3`, sex=`X2`, class=`X4`)
eps$X5 <- NULL
#
#tbl <- table(sub("^.*\\.", "", eps$ensembl_id.ver))
#writeLines(sprintf("EnsembleID version %s: N = %d", names(tbl),tbl))
for (j in 5:ncol(eps)) {
  eps[[j]] <- round(as.numeric(eps[[j]]), digits=3)
}
eps$gene <- sub("^;", "", eps$gene)
eps$gene <- sub(";.*$", "", eps$gene) #Keep 1st symb only
### Currently cannot handle ambiguous gene symbols (Switch to ENSG IDs!).
eps <- eps[!duplicated(eps[,c("gene","sex")]),] #REMOVES AMBIGUOUS GSYMs; TO BE ADDRESSED.
#
#F and M gene sets should be same.
eps_f <- eps[eps$sex=="F",]
eps_m <- eps[eps$sex=="M",]
if (!setequal(eps_f$gene, eps_m$gene)) {
  eps_f <- eps_f[eps_f$gene %in% intersect(eps_f$gene,eps_m$gene),]
  eps_m <- eps_m[eps_m$gene %in% intersect(eps_m$gene,eps_m$gene),]
  eps <- eps[eps$gene %in% intersect(eps_f$gene,eps_m$gene),]
}
eps <- eps[order(eps$gene, eps$sex),]
eps_f <- eps_f[order(eps_f$gene),]
eps_m <- eps_m[order(eps_m$gene),]
#
###
# From http://uswest.ensembl.org/biomart/
###
ensg2ncbi <- read_delim("data/biomart_ENSG2NCBI.tsv.gz", "\t", col_names = c("ensg","ensg.ver","ncbi","hgnc"), skip=1)
ensg2ncbi <- ensg2ncbi[!is.na(ensg2ncbi$ncbi),c(2,3,4)]
eps <- merge(eps, ensg2ncbi, by="ensg.ver", all.x=T, all.y=F)
writeLines(sprintf("Genes mapped to NCBI gene IDs: %d ; unmapped: %d ; pct: %.1f%%", 
                   nrow(eps[!is.na(eps$ncbi),]), nrow(eps[is.na(eps$ncbi),]), 
                   100*nrow(eps[!is.na(eps$ncbi),])/nrow(eps)))
write_csv(eps, path = "data/exfiles_eps.csv")
###
tissues <- read_csv("data/exfiles_tissues.csv")
cat("TISSUES:")
writeLines(sprintf("%d. %s", tissues$tissue_id, tissues$tissue_name))
#
### Gene-gene correlations from Gio.
ggc <- read_delim("data/result_recall.csv.gz", "\t", escape_double=F, trim_ws=T)
ggc <- rename(ggc, GiClass = `Gi class`, GjClass = `Gj class`, ABC = `dAUC`)
ggc$Gi <- sub("^;", "", ggc$Gi)
ggc$Gj <- sub("^;", "", ggc$Gj)
ggc$Gi <- sub(";.*$", "", ggc$Gi) #Keep 1st symb only
ggc$Gj <- sub(";.*$", "", ggc$Gj) #Keep 1st symb only


gi <- unique(ggc[,c("Gi","GiClass")])
gi <- rename(gi, sym = Gi, class = GiClass)
gi <- gi[order(gi$sym),]
gj <- unique(ggc[,c("Gj","GjClass")])
gj <- rename(gj, sym = Gj, class = GjClass)
gj <- gj[order(gj$sym),]


gene <- rbind(gi, gj)
gene <- unique(gene)
gene <- gene[order(gene$sym),]
#
tcrd <- read_csv("~/projects/idg/TCRD/data/tcrd_targets.csv")
tcrd <- rename(tcrd, sym=`protein.sym`, 
  uniprot=`protein.uniprot`, up_version=`protein.up_version`, ensemblid=`protein.stringid`, geneid=`protein.geneid`,
	name=`target.name`, fam=`target.fam`, tdl=`target.tdl`, idg2=`target.idg2`)
#
gene <- merge(gene, tcrd[,c("sym", "geneid", "uniprot", "up_version", "ensemblid", "name", "fam",  "tdl", "idg2")], 
           by="sym", all.x=T, all.y=F)
#
writeLines(sprintf("Genes mapped to TCRD: %d ; unmapped: %d (%.1f%%)", 
                   nrow(gene[!is.na(gene$geneid),]),
                   nrow(gene[is.na(gene$geneid),]),
                   100*nrow(gene[!is.na(gene$geneid),])/nrow(gene)))
#
writeLines(sprintf("GENE COUNT: %d", nrow(gene)))
write_csv(gene, path = "data/exfiles_gene.csv")


ggc$GiClass <- NULL #Save space
ggc$GjClass <- NULL #Save space
ggc$rho <- round(ggc$rho, digits=2)
ggc$RMSD <- round(ggc$RMSD, digits=2)
ggc$ABC <- round(ggc$ABC, digits=2)
ggc <- ggc[,c(1,2,4,3,5,6)] #reorder cols, metrics last: Gi,Gj,Cluster,rho,RMSD,ABC
ggc$wrho <- NA
ggc$tmoto <- NA
#
###
# Populate matrix with each named row an expression profile.  M_F mean.
eps_mx <- (as.matrix(eps_f[,5:44]) + as.matrix(eps_m[,5:44]))/2 #M+F mean
rownames(eps_mx) <- eps_f$gene
###
n_chunk <- 1e5
t0 <- proc.time()
i <- 1
while (i<nrow(ggc)) {
  i_next <- min(i+n_chunk, nrow(ggc)+1)
  gAs <- ggc$Gi[i:(i_next-1)]
  gBs <- ggc$Gj[i:(i_next-1)]
  epAs <- eps_mx[gAs,]
  epBs <- eps_mx[gBs,]

  ggc$tmoto[i:(i_next-1)] <- Tanimoto_mx(epAs, epBs)
  ggc$wrho[i:(i_next-1)] <- wPearson_mx(epAs, epBs)

  if ((i %% n_chunk)==1) {
    writeLines(sprintf("Progress: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc), 100*i/nrow(ggc), hf_utils$NiceTime((proc.time()-t0)[3])))
  }
  i <- i_next
}
writeLines(sprintf("Total: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc), 100*i/nrow(ggc), hf_utils$NiceTime((proc.time()-t0)[3])))
###
ggc$tmoto <- round(ggc$tmoto, digits=2)
ggc$wrho <- round(ggc$wrho, digits=2)
#
#Cluster "H+" and "H-" means rho>.7 and rho<-.7, so is redundant.
ggc$Cluster <- sub("H[+-] \\((.*)\\)", "\\1", ggc$Cluster)
tbl <- table(ggc$Cluster)
writeLines(sprintf("%12s: %6d", names(tbl), tbl))
#
tbl <- table(gene$target.fam)
writeLines(sprintf("%12s: %6d", names(tbl), tbl))
###
# TO LIMIT MEMORY FOR SHINYAPPS.IO (1G for Starter) MAY NEED TO REMOVE COLS HERE.
###
gzout <- gzfile("data/exfiles_ggc.csv.gz","w")
write.csv(ggc, gzout, row.names=F)
close(gzout)
###

###
# (1) Test with exfiles_plot_test.R.
# (2) Copy data files to Shiny app exfiles directory.
###
#
