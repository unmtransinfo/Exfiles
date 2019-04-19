#!/usr/bin/env Rscript
###
### Compute SABV measure[s] for input gene list.
###
library(readr)
library(data.table)
library(wCorr)
library(plotly)
#
Ruzicka <- function(A,B) {
  sum(pmin(A,B), rm.na=T)/sum(pmax(A,B), rm.na=T)
}
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  wCorr::weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
indata <- read_delim("data/Copy_of_RelevantMvsF_ADRs.tsv", "\t")
setDT(indata)
###
# Compute SABV measures for single genes, F vs. M.
###
eps <- read_delim("data/exfiles_eps.tsv", "\t", col_types=cols(.default=col_number(), ENSG=col_character(), SEX=col_character()))
setDT(eps)
#
###
#
writeLines(sprintf("Input ENSGs mapped to Exfiles: %d / %d", length(intersect(indata$ensg, eps$ENSG)), uniqueN(indata$ensg)))
#

### Expression profile for gene, sex.
expro <- function(ensg_this, sex_this, eps) {
  if (nrow(eps[ENSG==ensg_this & SEX==sex_this])==0) {
    message(sprintf("DEBUG: %s (%s) not found.", ensg_this, sex_this))
    return(NA)
  }
  return(as.numeric(eps[ENSG==ensg_this & SEX==sex_this, 3:ncol(eps)]))
}

n_err <- 0
for (ensg_this in unique(indata$ensg)) {
  expro_f <- expro(ensg_this, "F", eps)
  expro_m <- expro(ensg_this, "M", eps)
  if (is.na(expro_f) | is.na(expro_m)) {
    n_err <- n_err + 1
    next 
  }
  ruz <- Ruzicka(expro_f, expro_m)
  rho <- wPearson(expro_f, expro_m)
  #message(sprintf("DEBUG: ensg=%s ; ruz=%.3f ; rho=%.3f", ensg_this, ruz, rho))
  indata[ensg==ensg_this, exfiles_ruzicka_f_vs_m := ruz]
  indata[ensg==ensg_this, exfiles_rho_f_vs_m := rho]
  
}
message(sprintf("n_err: %d", n_err))

plots <- list()
plots[["Ruzicka"]] <- plot_ly(indata, type="histogram", x=~exfiles_ruzicka_f_vs_m, name="Ruzicka")
plots[["wPearson"]] <- plot_ly(indata, type="histogram", x=~exfiles_rho_f_vs_m, name="wPearson")
subplot(plots, nrows=1) %>%
  layout(title="Exfiles similarity measures", margin=list(t=120))


qtl <- quantile(indata$exfiles_ruzicka_f_vs_m, probs=seq(0, 1, .1), na.rm=T)
message(sprintf("Ruzicka: %5sile: %.3f\n", names(qtl), qtl))
qtl <- quantile(indata$exfiles_rho_f_vs_m, probs=seq(0, 1, .1), na.rm=T)
message(sprintf("wRho: %5sile: %.3f\n", names(qtl), qtl))

indata[, exfiles_combo_f_vs_m := exfiles_ruzicka_f_vs_m * exfiles_rho_f_vs_m]
setorder(indata, exfiles_combo_f_vs_m, na.last=T)

write_delim(indata, "data/Copy_of_RelevantMvsF_ADRs_ExfilesSABV.tsv", delim="\t")
###

