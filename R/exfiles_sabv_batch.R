#!/usr/bin/env Rscript
###
### Compute SABV measure[s] for input gene list.
###
library(readr)
library(data.table)
#
Ruzicka <- function(A,B) {
  sum(pmin(A,B), rm.na=T)/sum(pmax(A,B), rm.na=T)
}
###
indata <- read_delim("data/Copy_of_RelevantMvsF_ADRs.tsv", "\t")
setDT(indata)
###
# Compute SABV measures for single genes, F vs. M.
###
eps <- read_delim("data/exfiles_eps.csv", ",", col_types = cols(.default = col_number(), ensg=col_character(),
                                                                ncbi=col_character(), gene=col_character(),
                                                                sex=col_character()))
setDT(eps)
#
###
#
writeLines(sprintf("Input ENSGs mapped to Exfiles: %d / %d", length(intersect(indata$ensg, eps$ensg)), uniqueN(indata$ensg)))
#
#hugo_ids <- read_delim("data/hugo_protein-coding_gene.tsv", "\t", col_types=cols(.default=col_character()))
#setDT(hugo_ids)
#hugo_ids <- unique(hugo_ids[, .(hugo_symbol = symbol, hugo_uniprot = uniprot_ids, hugo_ensg = ensembl_gene_id)])
#
#writeLines(sprintf("Input UNIPROTs mapped to HUGO: %d / %d", length(intersect(indata$uniprot, hugo_ids$hugo_uniprot)), uniqueN(indata$uniprot)))
#writeLines(sprintf("HUGO ENSGs mapped to Exfiles: %d / %d", length(intersect(hugo_ids$hugo_ensg, eps$ensg)), uniqueN(eps$ensg)))
#
#indata <- merge(indata, hugo_ids, by.x="uniprot", by.y="hugo_uniprot", all.x=T, all.y=F)
#
#writeLines(sprintf("HUGO ENSGs mapped to Input and Exfiles: %d / %d", length(intersect(indata$hugo_ensg, eps$ensg)), uniqueN(indata$hugo_ensg)))
#

### Expression profile for gene, sex.
expro <- function(ensg_this, sex_this, eps) {
  #if (nrow(eps[ensg==ensg_this & sex==sex_this])>1) {
  #  #writeLines(sprintf("DEBUG: %s (%s) ambiguous: %s", ensg_this, sex_this, paste(eps[ensg==ensg_this & sex==sex_this, ensg], collapse=",")))
  #  return(NA)
  #}
  if (nrow(eps[ensg==ensg_this & sex==sex_this])==0) {
    writeLines(sprintf("DEBUG: %s (%s) not found.", ensg_this, sex_this))
    return(NA)
  }
  return(as.numeric(eps[ensg==ensg_this & sex==sex_this, 5:ncol(eps)]))
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
  #message(sprintf("DEBUG: ensg=%s ; ruz=%f", ensg_this, ruz))
  indata[ensg==ensg_this, exfiles_ruzicka_f_vs_m := ruz]
}
message(sprintf("n_err: %d", n_err))
hist(indata$exfiles_ruzicka_f_vs_m)

qtl <- quantile(indata$exfiles_ruzicka_f_vs_m, probs=seq(0, 1, .1), na.rm=T)
message(sprintf("%5sile: %.3f\n", names(qtl), qtl))

setorder(indata, exfiles_ruzicka_f_vs_m, na.last=T)

write_delim(indata, "data/Copy_of_RelevantMvsF_ADRs_ExfilesSABV.tsv", delim="\t")
###

