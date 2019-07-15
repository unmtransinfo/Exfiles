#!/usr/bin/env Rscript
###
library(readr)
library(data.table)


###

message(sprintf("=== TCRD-GTEx export for comparison:"))
ifile_tcrd <- "data/tcrd_expression-gtex.tsv.gz"
message(sprintf("%s", ifile_tcrd))

tcrd_gtex <- read_delim(ifile_tcrd, "\t", col_types=cols(.default=col_character(), number_value=col_double()), na=c("", "NULL"))
setDT(tcrd_gtex)

message(sprintf("Total GTEx data from TCRD table expression: %d", nrow(tcrd_gtex)))

for  (tag in names(tcrd_gtex)) {
  message(sprintf("Non-NA values: %9d [%s]", sum(!is.na(tcrd_gtex[[tag]])), tag))
}

message(sprintf("Unique genes: %d", uniqueN(tcrd_gtex[, protein_id])))
message(sprintf("Unique tissues: %d", uniqueN(tcrd_gtex[, tissue])))

tbl <- table(tcrd_gtex$sex, useNA="ifany")
message(sprintf("%6s: %9d\n", names(tbl), as.integer(tbl)))

tbl <- table(tcrd_gtex$age, useNA="ifany")
message(sprintf("%6s: %9d\n", names(tbl), as.integer(tbl)))

tbl <- table(tcrd_gtex$qual_value, useNA="ifany")
message(sprintf("%6s: %9d\n", names(tbl), as.integer(tbl)))

qtl <- quantile(tcrd_gtex$number_value, seq(0, 1, .1), na.rm=T)
message(sprintf("%sile: %.1f\n", names(qtl), qtl))

hist(tcrd_gtex$number_value, main="GTEx expression histogram")

###

message(sprintf("=== Exfiles-GTEx SABV workflow output file:"))
ifile_sabv <- "data/gtex_rnaseq_sabv_alltissues.tsv.gz"
message(sprintf("%s", ifile_sabv))

gtex_sabv <- read_delim(ifile_sabv, "\t", col_types=cols(.default=col_double(), ENSG=col_character(), SEX=col_character(), SMTSD=col_character(), TPM_LEVEL=col_character(), TPM_LEVEL_BYSEX=col_character()), na=c("", "NULL"))
setDT(gtex_sabv)

message(sprintf("Total GTEx data from Exfiles workflow: %d", nrow(gtex_sabv)))
message(sprintf("Unique genes: %d", uniqueN(gtex_sabv[, ENSG])))
message(sprintf("Unique tissues: %d", uniqueN(gtex_sabv[, SMTSD])))

message(sprintf("Tissue in this dataset and NOT in TCRD file: \"%s\"\n", setdiff(gtex_sabv$SMTSD, tcrd_gtex$tissue)))


for (tag in c("SEX", "TPM_LEVEL", "TPM_LEVEL_BYSEX")) {
  tbl <- table(gtex_sabv[[tag]], useNA="ifany")
  message(sprintf("%s: %6s: %9d\n", tag, names(tbl), as.integer(tbl)))
}

