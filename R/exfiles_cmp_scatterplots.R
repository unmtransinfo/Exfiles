#!/usr/bin/env Rscript
###
###
library(readr)
library(dplyr)
library(plotly, quietly=T)


###
ggc <- read_delim("data/exfiles_ggc.tsv", "\t", col_types="cccdd")

## Assure Ga<Gb (alphabetical).  This allows merging on both columns.
ggc["Ga_tmp"] <- mapply(min, ggc$Ga, ggc$Gb, USE.NAMES=F)
ggc["Gb_tmp"] <- mapply(max, ggc$Ga, ggc$Gb, USE.NAMES=F)
ggc$Ga <- NULL
ggc$Gb <- NULL
ggc <- rename(ggc, Ga=Ga_tmp, Gb=Gb_tmp, rho=wRho, sim=Ruzicka)


for (group in c("F","M","FM"))
{
  qtl <- quantile(ggc$rho[ggc$Cluster==group], probs = c(seq(0,.09,.01), seq(0.1, 1, 0.1)))
  writeLines(sprintf("Rho_%s: %5s-ile: %.3f", group, names(qtl), qtl))
  cat("\n")
}

### Compare the metrics.
writeLines(sprintf("Ruzicka vs. WPearson: %.2f", cor(ggc$sim, ggc$rho)))


### Scatter plots
plots <- list()
plots[["Ruzicka vs. WPearson"]] <- plot_ly(type="scatter", mode="markers", 
                                           data=ggc[sample(nrow(ggc),5e3),], x=~rho, y=~sim) %>%
  layout(xaxis=list(title="wRho"), yaxis=list(title="Ruzicka"))
#
subplot(plots, nrows=1, margin=0.05, shareX=F, shareY=F, titleX=T, titleY=T) %>%
  layout(title="Ex-files: Ruzicka similarity vs. Weighted Pearson correlation", margin=list(t=80, b=60, l=30),
         font=list(family="Arial", size=14), showlegend=F)
#
