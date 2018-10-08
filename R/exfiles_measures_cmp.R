#!/usr/bin/env Rscript
###
#
library(readr)
library(dplyr)
library(wCorr)
library(plotly, quietly=T)
###
#
RMSD <- function(A,B) {
  sqrt(sum((A-B)^2)/length(A))
}
#
Euclid <- function(A,B) {
  sqrt(sum((A-B)^2))
}
#
Tanimoto <- function(A,B) { #Vector version (slow)
  (A%*%B) / (A%*%A + B%*%B - A%*%B)
}
#
wPearson <- function(A,B) { #Vector version (slow)
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
#
Ruzicka <- function(A,B) {
  sum(pmin(A,B))/sum(pmax(A,B))
}
###
eps <- read_delim("data/exfiles_eps.tsv", "\t", col_types=cols(SEX=col_character()))
eps <- eps[eps$SEX=="F",] #arbitrary
eps$SEX <- NULL
# Random sample:
N <- 200
eps <- eps[sample(nrow(eps),N),]
#
ggc <- data.frame(ENSGA = rep(eps$ENSG, each=N), ENSGB=rep(eps$ENSG, times=N))
ggc$ENSGA <- as.character(ggc$ENSGA)
ggc$ENSGB <- as.character(ggc$ENSGB)
ggc <- ggc[ggc$ENSGA<ggc$ENSGB,]
ggc <- ggc[order(ggc$ENSGA, ggc$ENSGB),]
rownames(ggc) <- NULL
ggc['wRho'] <- NA
ggc['Ruzicka'] <- NA
ggc['Tanimoto'] <- NA
ggc['RMSD'] <- NA
ggc['Euclid'] <- NA
#
###
# Compute
#
ensgA <- ""
for (i in 1:nrow(ggc)) {
  if (ggc$ENSGA[i]!=ensgA) {
    ensgA <- ggc$ENSGA[i]
    epsA <- as.numeric(eps[eps$ENSG==ensgA,2:ncol(eps)])
  }
  epsB <- as.numeric(eps[eps$ENSG==ggc$ENSGB[i],2:ncol(eps)])
  ggc$wRho[i] <- wPearson(epsA,epsB)
  ggc$Ruzicka[i] <- Ruzicka(epsA,epsB)
  ggc$RMSD[i] <- RMSD(epsA,epsB)
  ggc$Euclid[i] <- Euclid(epsA,epsB)
  ggc$Tanimoto[i] <- Tanimoto(epsA,epsB)
}
writeLines(sprintf("Computed: %d", i))


### Compare the metrics.
writeLines(sprintf("Ruzicka vs. WPearson: %.2f", cor(ggc$Ruzicka, ggc$wRho)))


### Scatter plots
plots <- list()
#
plots[["Ruzicka_vs_wRho"]] <- plot_ly(type="scatter", mode="markers", marker=list(symbol="dot", size=4),
	data=ggc, x=~wRho, y=~Ruzicka) %>%
	layout(yaxis=list(title="Ruzicka"))
plots[["Ruzicka_vs_Euclid"]] <- plot_ly(type="scatter", mode="markers", marker=list(symbol="dot", size=4),
	data=ggc, x=~Euclid, y=~Ruzicka) %>%
	layout(yaxis=list(title="Ruzicka"))
plots[["Ruzicka_vs_RMSD"]] <- plot_ly(type="scatter", mode="markers", marker=list(symbol="dot", size=4),
	data=ggc, x=~RMSD, y=~Ruzicka) %>%
	layout(yaxis=list(title="Ruzicka"))
plots[["Ruzicka_vs_Tanimoto"]] <- plot_ly(type="scatter", mode="markers", marker=list(symbol="dot", size=4),
	data=ggc, x=~Tanimoto, y=~Ruzicka) %>%
	layout(yaxis=list(title="Ruzicka"))
#
subplot(plots, nrows=2, margin=0.05, shareX=F, shareY=T, titleX=F, titleY=T) %>%
  layout(title="Ex-files: Measures comparisons", margin=list(t=80, b=60, l=30),
         font=list(family="Arial", size=14), showlegend=F) %>%
  add_annotations(names(plots), x=c(.1, .8, .1, .8), y=c(1, 1, .4, .4), xref="paper", yref="paper", showarrow=F)
#
