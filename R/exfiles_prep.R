#!/usr/bin/env Rscript
###
### This processes the original files provided by Giovanni (early Feb 2018): 
### (1) result_recall.csv: gene-gene associations
### (2) receptors_all.csv: expression profiles
###
library(readr)
library(wCorr)
library(dplyr)
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
# Count and order must correspond with tpm columns.
tissues <- read_delim("data/tissue_order.csv",  ";", col_names = c("tissue_name", "tissue_id"))
n_tissues <- nrow(tissues)
cat("TISSUES:")
writeLines(sprintf("%d. %s", tissues$tissue_id, tissues$tissue_name))
write_csv(tissues[,c("tissue_id","tissue_name")], path = "data/exfiles_tissues.csv")
#
### Expresion profiles from Gio.
### LOG10(1+TPM)
#eps <- read_delim("data/receptors_all.csv.gz", "\t", col_names=F)
#eps <- transform(eps, tissue=do.call(rbind, strsplit(X5, '[, ]+', fixed=F)), stringsAsFactors=F)
#eps <- rename(eps, ensg.ver=`X1`, gene=`X3`, sex=`X2`, class=`X4`)
#eps$X5 <- NULL
#
eps <- read_delim("data/tpm_vectors.csv", "\t", col_names=F)
eps <- transform(eps, tissue=do.call(rbind, strsplit(X3, ';', fixed=F)), stringsAsFactors=F)
eps <- dplyr::rename(eps, ensg=`X1`, sex=`X2`)
eps$X3 <- NULL
#
if (ncol(eps)-2  != n_tissues) {
  stop(sprintf("TPM column count vs. tissue count mismatch (%d != %d)", ncol(eps)-2, n_tissues))
}
#
#tbl <- table(sub("^.*\\.", "", eps$ensembl_id.ver))
#writeLines(sprintf("EnsembleID version %s: N = %d", names(tbl),tbl))
for (j in 3:ncol(eps)) {
  eps[[j]] <- round(as.numeric(eps[[j]]), digits=2)
}
#eps$gene <- sub("^;", "", eps$gene)
#eps$gene <- sub(";.*$", "", eps$gene) #Keep 1st symb only
### Currently cannot handle ambiguous gene symbols (Switch to ENSG IDs!).
#eps <- eps[!duplicated(eps[,c("gene","sex")]),] #REMOVES AMBIGUOUS GSYMs; TO BE ADDRESSED.
#

###
# From http://uswest.ensembl.org/biomart/
###
ensg2ncbi <- read_delim("data/biomart_ENSG2NCBI.tsv", "\t", col_names = c("ensg","ensg.ver","ncbi","hgnc"), skip=1)
ensg2ncbi <- ensg2ncbi[!is.na(ensg2ncbi$ncbi),c("ensg","ncbi","hgnc")]
eps <- merge(eps, ensg2ncbi, by="ensg", all.x=T, all.y=F)
writeLines(sprintf("Genes mapped to NCBI gene IDs: %d ; unmapped: %d ; pct: %.1f%%", 
                   nrow(eps[!is.na(eps$ncbi),]), nrow(eps[is.na(eps$ncbi),]), 
                   100*nrow(eps[!is.na(eps$ncbi),])/nrow(eps)))
eps <- dplyr::rename(eps, gene = hgnc)
eps <- eps[,c("ensg", "ncbi", "gene", "sex", paste0("tissue.", 1:n_tissues))]
colnames(eps) <- c("ensg", "ncbi", "gene", "sex", tissues$tissue_name)
eps <- eps[!is.na(eps$gene),]
#
#F and M gene sets should be same.
eps_f <- eps[eps$sex=="F",]
eps_m <- eps[eps$sex=="M",]
eps_f$sex <- NULL
eps_m$sex <- NULL
#
#if (!setequal(eps_f$ensg, eps_m$ensg)) {
#  eps_f <- eps_f[eps_f$ensg %in% intersect(eps_f$ensg,eps_m$ensg),]
#  eps_m <- eps_m[eps_m$ensg %in% intersect(eps_m$ensg,eps_m$ensg),]
#  eps <- eps[eps$ensg %in% intersect(eps_f$ensg,eps_m$ensg),]
#}
#eps <- eps[order(eps$ensg, eps$sex),]
#eps_f <- eps_f[order(eps_f$ensg),]
#eps_m <- eps_m[order(eps_m$ensg),]
#
# Unfortunately must index by gene symbols for now.
if (!setequal(eps_f$gene, eps_m$gene)) {
  eps_f <- eps_f[eps_f$gene %in% intersect(eps_f$gene,eps_m$gene),]
  eps_m <- eps_m[eps_m$gene %in% intersect(eps_m$gene,eps_m$gene),]
  eps <- eps[eps$gene %in% intersect(eps_f$gene,eps_m$gene),]
}
eps <- eps[!duplicated(eps[,c("gene","sex")]),]
eps_f <- eps_f[!duplicated(eps_f$gene),]
eps_m <- eps_m[!duplicated(eps_m$gene),]
#
write_csv(eps, path = "data/exfiles_eps.csv")
#
### Gene-gene correlations from Gio.
# ID-ed by gene symbols (sub-optimal).
ggc <- read_delim("data/result_recall.csv.gz", "\t", escape_double=F, trim_ws=T)
ggc <- rename(ggc, GiClass = `Gi class`, GjClass = `Gj class`, ABC = `dAUC`)
ggc$Gi <- sub("^;", "", ggc$Gi)
ggc$Gj <- sub("^;", "", ggc$Gj)
ggc$Gi <- sub(";.*$", "", ggc$Gi) #Keep 1st symb only
ggc$Gj <- sub(";.*$", "", ggc$Gj) #Keep 1st symb only
ggc$GiClass <- NULL #Save space
ggc$GjClass <- NULL #Save space

gi <- unique(ggc[,c("Gi")])
gi <- rename(gi, sym = Gi)
gi <- gi[order(gi$sym),]
gj <- unique(ggc[,c("Gj")])
gj <- rename(gj, sym = Gj)
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
# May be multiple rows per gene symbol. (Switch IDs to ENSG!)
gene <- gene[!duplicated(gene$sym),]
#
writeLines(sprintf("Genes mapped to TCRD: %d ; unmapped: %d (%.1f%%)", 
                   nrow(gene[!is.na(gene$geneid),]),
                   nrow(gene[is.na(gene$geneid),]),
                   100*nrow(gene[!is.na(gene$geneid),])/nrow(gene)))
#
writeLines(sprintf("GENE COUNT: %d", nrow(gene)))
write_csv(gene, path = "data/exfiles_gene.csv")
#
#
ggc$rho <- round(ggc$rho, digits=2)
ggc$RMSD <- round(ggc$RMSD, digits=2)
ggc$ABC <- round(ggc$ABC, digits=2)
ggc <- ggc[,c(1,2,4,3,5,6)] #reorder cols, metrics last: Gi,Gj,Cluster,rho,RMSD,ABC
ggc[["wrho"]] <- NA
#ggc[["tmoto"]] <- NA
#
###
# Remove rows from ggc for which profiles do not exist (due to symbol ambiguity).
ggc <- ggc[ggc$Gi %in% eps$gene & ggc$Gj %in% eps$gene,]
###
# Populate matrix with each named row an expression profile.  M_F mean.
eps_mx <- (as.matrix(eps_f[,4:ncol(eps_f)]) + as.matrix(eps_m[,4:ncol(eps_m)]))/2 #M+F mean
rownames(eps_mx) <- eps_f$gene
#
# Must have profiles for each gene in ggc.
writeLines(sprintf("DEBUG: Expression profiles (F): %d", length(unique(eps_f$gene))))
writeLines(sprintf("DEBUG: Expression profiles (M): %d", length(unique(eps_m$gene))))
writeLines(sprintf("DEBUG: Expression profiles matrix ((F+M)/2): %d", nrow(eps_mx)))
writeLines(sprintf("DEBUG: ggc Gi genes: %d", length(unique(ggc$Gi))))
writeLines(sprintf("DEBUG: ggc Gj genes: %d", length(unique(ggc$Gj))))
###
#for (g in c(ggc$Gi,ggc$Gj)) {
#  if (!(g %in% rownames(eps_mx))) {
#    writeLines(sprintf("DEBUG: %s not in eps_mx", g))
#    break
#  }
#}
###
source(paste0(Sys.getenv("HOME"),"/lib/R/time_utils.R"))
n_chunk <- 1e5
t0 <- proc.time()
i <- 1
while (i<nrow(ggc)) {
  i_next <- min(i+n_chunk, nrow(ggc)+1)
  gAs <- ggc$Gi[i:(i_next-1)]
  gBs <- ggc$Gj[i:(i_next-1)]
  epAs <- eps_mx[gAs,]
  epBs <- eps_mx[gBs,]

  #ggc$tmoto[i:(i_next-1)] <- Tanimoto_mx(epAs, epBs)
  ggc$wrho[i:(i_next-1)] <- wPearson_mx(epAs, epBs)

  if ((i %% n_chunk)==1) {
    writeLines(sprintf("Progress: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc), 100*i/nrow(ggc), time_utils$NiceTime((proc.time()-t0)[3])))
  }
  i <- i_next
}
writeLines(sprintf("Total: %7d / %7d (%.1f%%) ; elapsed: %s", i-1, nrow(ggc), 100*i/nrow(ggc), time_utils$NiceTime((proc.time()-t0)[3])))
###
#ggc$tmoto <- round(ggc$tmoto, digits=2)
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
