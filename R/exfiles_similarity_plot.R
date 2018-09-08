#!/usr/bin/env Rscript
#############################################################################
library(readr)
library(plotly)

exfiles_pearson <- read_delim("data/gtex_rnaseq_prep_profiles_PearsonNxN.tsv", "\t")

exfiles_spearman <- read_delim("data/gtex_rnaseq_prep_profiles_SpearmanNxN.tsv", "\t")

exfiles_tanimoto <- read_delim("data/gtex_rnaseq_prep_profiles_TanimotoNxN.tsv", "\t")

exfiles <- merge(exfiles_pearson, exfiles_spearman, by=c("ENSGA","SEXA","ENSGB","SEXB"), all=T)

exfiles <- merge(exfiles, exfiles_tanimoto, by=c("ENSGA","SEXA","ENSGB","SEXB"), all=T)

###
for (field in c("Pearson", "SpearmanRho", "Tanimoto")) {
  print(sprintf("%s mean: %f median: %f", field, mean(exfiles[[field]], na.rm=T), median(exfiles[[field]], na.rm=T)))
  qs <- quantile(exfiles[[field]], probs = c(seq(.1, .9, .1), seq(.91, .99, .01), seq(.991, 1, .001)), na.rm=T)
  for (name in names(qs))
  {
    print(sprintf("%s %s quantile: %f", field, name, qs[[name]]))
  }
}


###
plot_ly(exfiles[sample(nrow(exfiles), 1e5),], x=~Pearson, y=~SpearmanRho, type="scatter", mode="markers", marker = list(symbol="circle", size=2))

###
plots = list()
plots[["Pearson"]] = plot_ly(exfiles, y=~Pearson, type="box",
                             name="Pearson",
                             boxmean=T, boxpoints=F,
                             marker=list(symbol="dot", opacity=1))
plots[["Spearman"]] = plot_ly(exfiles, y=~SpearmanRho, type="box",
                              name="Spearman",
                             boxmean=T, boxpoints=F,
                             marker=list(symbol="dot", opacity=1))
plots[["Tanimoto"]] = plot_ly(exfiles, y=~Tanimoto, type="box",
                              name="Tanimoto",
                             boxmean=T, boxpoints=F,
                             marker=list(symbol="dot", opacity=1))
subplot(plots, nrows=1) %>%
  layout(title = paste0("GTEx expression profile comparison metrics"),
         margin = list(t=100),
         showlegend=F, font = list(family="Arial", size=14))

