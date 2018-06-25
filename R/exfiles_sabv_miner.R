#!/usr/bin/env Rscript
###
### Find genes with max SABV.
###
library(readr)
library(dplyr)
library(wCorr)
#library(diptest)
library(plotly, quietly=T)

RMSD <- function(A,B) {
  sqrt(sum((A-B)^2)/length(A))
}
Tanimoto <- function(A,B) {
  (A %*% B) / (A %*% A + B %*% B - A %*% B)
}
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
ABC <- function(A,B) {
  abc <- 0
  for (i in 1:(length(A)-1)) {
    Amid <- mean(A[i],A[i+1])
    Bmid <- mean(B[i],B[i+1])
    if (A[i]>=B[i]) {
      if (A[i+1]>=B[i+1]) {
        abc <- abc + (AreaUnderLineSegment(A[i], A[i+1], 1) - AreaUnderLineSegment(B[i], B[i+1], 1))
      } else {
        abc <- abc + (AreaUnderLineSegment(A[i], Amid, .5) - AreaUnderLineSegment(B[i], Bmid, .5))
        abc <- abc + (AreaUnderLineSegment(Bmid, B[i+1], .5) - AreaUnderLineSegment(Amid, A[i+1], .5))
      }
    } else { 
      if (A[i+1]<B[i+1]) {
        abc <- abc + (AreaUnderLineSegment(B[i], B[i+1], 1) - AreaUnderLineSegment(A[i], A[i+1], 1))
      } else {
        abc <- abc + (AreaUnderLineSegment(B[i], Bmid, .5) - AreaUnderLineSegment(A[i], Amid, .5))
        abc <- abc + (AreaUnderLineSegment(Amid, A[i+1], .5) - AreaUnderLineSegment(Bmid, B[i+1], .5))
      }
    }
  }
  return(abc)
}
#
AreaUnderLineSegment <- function(y1, y2, w) {
  a <- min(y1,y2) * w
  a <- a + 0.5 * w * abs(y1-y2)
  return(a)
}
###
# 40 GTEx tissues:
tissues <- read_delim("data/exfiles_tissues.csv", ",")
ntissue <- nrow(tissues)
writeLines(sprintf("Tissues: %d", ntissue))

###
ggc <- read_delim("data/exfiles_ggc.csv.gz", ",", col_types="cccddddd")

## Assure Gi<Gj (alphabetical).  This allows merging on both columns.
ggc$Gi_tmp <- mapply(min, ggc$Gi, ggc$Gj, USE.NAMES=F)
ggc$Gj_tmp <- mapply(max, ggc$Gi, ggc$Gj, USE.NAMES=F)
ggc$Gi <- NULL
ggc$Gj <- NULL
ggc <- rename(ggc, Gi=Gi_tmp, Gj=Gj_tmp)


for (group in c("F","M","MF"))
{
  qtl <- quantile(ggc$rho[ggc$Cluster==group], probs = c(seq(0,.09,.01), seq(0.1, 1, 0.1)))
  writeLines(sprintf("Rho_%s: %5s-ile: %.3f", group, names(qtl), qtl))
  cat("\n")
}

### Find gene-pairs where Rho-F exists but not Rho-M and not Rho-MF.
ggc_f <- rename(ggc[ggc$Cluster=="F",], rho_f=rho, RMSD_f=RMSD, ABC_f=ABC)
ggc_m <- rename(ggc[ggc$Cluster=="M",], rho_m=rho, RMSD_m=RMSD, ABC_m=ABC)
ggc_mf <- rename(ggc[ggc$Cluster=="MF",], rho_mf=rho, RMSD_mf=RMSD, ABC_mf=ABC)
ggc_f$Cluster <- NULL
ggc_m$Cluster <- NULL
ggc_mf$Cluster <- NULL
ggc_f_only <- merge(ggc_f,ggc_m, by=c("Gi","Gj"), all.x=T, all.y=F)
ggc_f_only <- merge(ggc_f_only,ggc_mf, by=c("Gi","Gj"), all.x=T, all.y=F)
ggc_f_only <- ggc_f_only[is.na(ggc_f_only$rho_m) | is.na(ggc_f_only$rho_mf),]
ggc_f_only <- ggc_f_only[order(-ggc_f_only$rho_f),]
for (colname in c("rho_m",   "RMSD_m", "ABC_m", "rho_mf", "RMSD_mf", "ABC_mf")) {
  ggc_f_only[[colname]] <- NULL
}


###
# Find single genes which have large F vs. M RMSD.
###
# Ambiguity issues: Symbol AARSD1 maps to ENSG00000108825.13 and ENSG00000266967.2 with different profiles.
###
eps <- read_csv("data/exfiles_eps.csv")
eps$ncbi <- NULL
eps$hgnc <- NULL
#eps$tissue.total <- rowSums(eps[,(ncol(eps)-ntissue+1):ncol(eps)])
#qtl <- quantile(eps$tissue.total, probs = c(seq(0,.09,.01), seq(0.1, 1, 0.1)))
#writeLines(sprintf("Total_expression: %5s-ile: %.3f", names(qtl), qtl))
#eps <- eps[eps$tissue.total>10,]  #Filter low total expression to avoid noise driven correlations.
eps <- eps[!duplicated(eps[,c("gene","sex")]),] #REMOVES AMBIGUOUS GSYMs; TO BE ADDRESSED.

eps_f <- eps[eps$sex=="F",]
eps_m <- eps[eps$sex=="M",]
#
###
# Individual tissue expression distributions.  With unimodality diptest.
for (i in 1:ntissue) {
  etissue_this <- eps[[sprintf("tissue.%d",i)]]
  writeLines(sprintf("%44s: Min: %4.2f ; Median: %4.2f ; Mean: %4.2f ; Max: %4.2f",
                     tissues$tissue_name[i], min(etissue_this), median(etissue_this), 
                     mean(etissue_this), max(etissue_this)))
}
###
#Expression histograms for each tissue:
xax=list(type="linear", title="")
yax=list(type="linear", title="")
plots <- list()
for (i in 1:ntissue) {
  plots[[tissues$tissue_name[i]]] <- plot_ly(type="histogram", x=eps[[sprintf("tissue.%d",i)]], 
        autobinx=F, xbins=list(start=0, end=3, size=.2), text=tissues$tissue_name[i])
}
#
n_rows <- 5
n_cols <- ceiling(as.numeric(length(plots))/n_rows)
n_cells <- n_rows*n_cols
xmx <- matrix(rep(seq(0, (n_cols-1)/n_cols, 1/n_cols), n_rows), byrow=T, ncol=n_cols)*1.1 #Hackfactor
ymx <- (-matrix(rep(seq(0, (n_rows-1)/n_rows, 1/n_rows), n_cols), byrow=F, ncol=n_cols)*1.1+1) #Hackfactor
subplot(plots, nrows=n_rows, margin=0.02, shareX=T, shareY=T, titleX=F, titleY=F) %>%
  layout(title="GTEx expression distributions by tissue", margin=list(t=80, b=60, l=30),
         xaxis=xax, yaxis=yax,
         font=list(family="Arial", size=14), showlegend=F) %>%
  add_annotations(text=paste0(1:n_cells, ". ", names(plots)),
    x=as.vector(t(xmx)), y=as.vector(t(ymx)),
	  xref="paper", yref="paper", align="left", font=list(family="Arial", size=14), showarrow=F) %>%
add_annotations(text="<I>UNITS: LOG<SUB>10</SUB>(1+TPM)</I>", x=.5, y=-.1, xref="paper", yref="paper", align="center", font=list(family="Arial", size=14), showarrow=F)
###
###
#Expression boxplots for each tissue:
epsx <- tidyr::gather(eps, key="tissue", value="expression", (ncol(eps)-ntissue+1):ncol(eps))
for (i in 1:ntissue) {
  epsx$tissue[epsx$tissue==sprintf("tissue.%d",i)] <- tissues$tissue_name[i]
}
#
plot_ly(epsx, type="box", x=~tissue, y=~expression, color=~sex, colors=c("#ff7777","#7777ff"),
        boxmean=T, boxpoints=F,  marker=list(symbol="dot",opacity=1)) %>%
  add_annotations(x=c(0.5), y=c(0.9),
                  xanchor="center",xref="paper",yref="paper",showarrow=F,
                  text=c(sprintf(""))) %>%
  layout(title="GTEx expression distributions by tissue", boxmode="group", 
         xaxis=list(title=""),yaxis=list(title="Expression (TPM)"), 
         margin=list(t=80,l=60,b=140,r=40),showlegend=T, legend = list(x=.9, y=1),
         font=list(family="Arial",size=12), titlefont=list(size=24))
rm(epsx)

###
genes <- intersect(eps_f$gene, eps_m$gene)
genes <- intersect(genes, union(ggc$Gi,ggc$Gj))

g_FvsM <- data.frame(gene = sort(genes), rho_fm=NA, wrho_fm=NA, rmsd_fm=NA, abc_fm=NA, tmoto_fm=NA)

### Expression profile for gene, sex.
expro <- function(gene,sex,eps) {
  eps_this <- eps[eps$gene==gene & eps$sex==sex,]
  if (nrow(eps_this)>1) {
    writeLines(sprintf("DEBUG: %s (%s) ambiguous: %s", gene, sex, paste(eps_this$ensg.ver, collapse=",")))
    return(NA)
  }
  return(as.numeric(eps_this[,5:44]))
}

for (i in 1:nrow(g_FvsM)) {
  gene <- g_FvsM$gene[i]
  expro_f <- expro(gene,"F",eps)
  expro_m <- expro(gene,"M",eps)
  if (is.na(expro_f) | is.na(expro_m)) { next }
  g_FvsM$rho_fm[i] <- cor(expro_f, expro_m, method="pearson")
  g_FvsM$rmsd_fm[i] <- RMSD(expro_f, expro_m)
  g_FvsM$abc_fm[i] <- ABC(expro_f, expro_m)
  g_FvsM$tmoto_fm[i] <- Tanimoto(expro_f, expro_m)
  g_FvsM$wrho_fm[i] <- wPearson(expro_f, expro_m)
}


###
### Compare the metrics.

writeLines(sprintf("ABC vs. RMSD: %.2f", cor(g_FvsM$abc_fm, g_FvsM$rmsd_fm)))
writeLines(sprintf("ABC vs. Tanimoto: %.2f", cor(g_FvsM$abc_fm, g_FvsM$tmoto_fm)))
writeLines(sprintf("Tanimoto vs. RMSD: %.2f", cor(g_FvsM$tmoto_fm, g_FvsM$rmsd_fm)))
writeLines(sprintf("Rho vs. RMSD: %.2f", cor(g_FvsM$rho_fm, g_FvsM$rmsd_fm, use="pairwise.complete.obs")))
writeLines(sprintf("Rho vs. weightedRho: %.2f", cor(g_FvsM$rho_fm, g_FvsM$wrho_fm, use="pairwise.complete.obs")))


### Scatter plots
plots <- list()
plots[["ABC vs. RMSD"]] <- plot_ly(type="scatter", mode="markers", data=g_FvsM, x=~rmsd_fm, y=~abc_fm, text=~gene)
plots[["ABC vs. Tanimoto"]] <- plot_ly(type="scatter", mode="markers", data=g_FvsM, x=~tmoto_fm, y=~abc_fm, text=~gene)
plots[["Tanimoto vs. RMSD"]] <- plot_ly(type="scatter", mode="markers", data=g_FvsM, x=~rmsd_fm, y=~tmoto_fm, text=~gene)
plots[["Rho vs. RMSD"]] <- plot_ly(type="scatter", mode="markers", data=g_FvsM, x=~rmsd_fm, y=~rho_fm, text=~gene)
#
subplot(plots, nrows=2, margin=0.05, shareX=F, shareY=F, titleX=F, titleY=F) %>%
  layout(title="Ex-files: Metrics", margin=list(t=80, b=60, l=30),
         font=list(family="Arial", size=14), showlegend=F) %>%
  add_annotations(text=names(plots), x=c(.2,.8,.2,.8), y=c(.9,.9,.4,.4), xref="paper", yref="paper", align="center", font=list(family="Arial", size=14), showarrow=F) %>%
add_annotations(text="<I>(All metrics same gene F,M)</I>", x=.5, y=-.1, xref="paper", yref="paper", align="center", font=list(family="Arial", size=14), showarrow=F)
#
plot_ly(type="scatter", mode="markers", data=g_FvsM, x=~rho_fm, y=~wrho_fm, text=~gene) %>%
  layout(title="Ex-files: Metrics: Weighted Correlation vs. Correlation", margin=list(t=80, b=60, l=30),
         font=list(family="Arial", size=14), showlegend=F)
#
#stop() #DEBUG

## What is sexy (max-SABV)?  Vive la difference!

#g_FvsM <- g_FvsM[order(g_FvsM$rho/g_FvsM$rmsd_fm),] 

g_FvsM <- g_FvsM[order(g_FvsM$tmoto_fm),]

#g_FvsM <- g_FvsM[order(-g_FvsM$abc_fm),]

## Plot SABV-max genes.


xaxis = list(title="", type="category", range=tissues$tissue_name, showticklabels=F)
yaxis = list(title="Expression", range=c(0,6))
#
n_rows <- 6
n_cols <- 4
plots <- list()
for (i in 1:(n_rows*n_cols)) {
  gene <- g_FvsM$gene[i]
  expro_f <- expro(gene,"F", eps)
  expro_m <- expro(gene,"M", eps)
  #
  rho <- cor(expro_f, expro_m, method="spearman", use="pairwise.complete.obs")
  rmsd <- RMSD(expro_f, expro_m)
  abc <- ABC(expro_f, expro_m)
  tmoto <- Tanimoto(expro_f, expro_m)
  #
  xref <- ifelse(i%%n_cols==1, "x", sprintf("x%d",(i-1)%%n_cols+1))
  yref <- ifelse(floor((i-1)/n_cols)==0, "y", sprintf("y%d", floor((i-1)/n_cols)+1))
  #
  plots[[i]] <- plot_ly(name=gene) %>%
    add_trace(name="F", x=tissues$tissue_name, y=expro_f, marker=list(color="red", size=10), line=list(color="red"),
        type="scatter", mode="lines+markers", legendgroup="1") %>%
    add_trace(name="M", x=tissues$tissue_name, y=expro_m, marker=list(color="blue", size=10), line=list(color="blue"),
              type="scatter", mode="lines+markers", legendgroup="1") %>%
  layout(xaxis=xaxis, yaxis=yaxis,
         annotations=list(text=sprintf("%s<br>Rho: %.2f ; Tanimoto: %.2f<br>ABC: %.2f ; RMSD: %.2f", gene, rho, tmoto, abc, rmsd), 
           x=ntissue, y=14, align="center", xref=xref, yref=yref, showarrow=F), 
           font=list(family="Arial", size=12))
  #
}
#
p <- subplot(plots, nrows=n_rows, margin=0.01, shareX=T, shareY=T, titleX=F, titleY=F) %>%
  layout(title="Ex-files: SABV Mining<br>", margin=list(t=80, b=80, l=30),
  font=list(family="Arial", size=14), showlegend=F) %>%
  add_annotations(text="Expression vs. Tissues (F:red ; M:blue) <I>Vive la difference!</I>", x=.5, y=-.1, xref="paper", yref="paper", align="center", font=list(family="Arial", size=16), showarrow=F)
#
p
#
