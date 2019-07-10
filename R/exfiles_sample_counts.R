#!/usr/bin/env Rscript
###
library(readr)
library(data.table)
#
# Columns: ENSG    SMTSD   SAMPID  SMATSSCR        SEX     AGE     DTHHRDY TPM
fpath <- "data/gtex_rnaseq_prep_sample.tsv"
samp <- read_delim(fpath, "\t")
setDT(samp)

writeLines(sprintf("TPM FILE: %s", fpath))
writeLines(sprintf("Columns: %s", paste0(colnames(samp), collapse=", ")))
writeLines(sprintf("Total TPMs: %d", nrow(samp)))
writeLines(sprintf("Total unique gene IDs: %d", length(unique(samp$ENSG))))
tbl <- table(samp$SEX)
writeLines(sprintf("TPMs (SEX = %s): %d (%.1f%%)", names(tbl), tbl, 100*tbl/nrow(samp)))
#
samp_samp <- unique(samp[,.(SAMPID,SEX)])
writeLines(sprintf("Total unique sample IDs: %d", length(unique(samp_samp$SAMPID))))
tbl <- table(samp_samp$SEX)
writeLines(sprintf("Samples (SEX = %s): %d (%.1f%%)", names(tbl), tbl, 100*tbl/nrow(samp_samp)))
#
samp[["SUBJID"]] <- sub("^([^-]+-[^-]+)-.*$", "\\1", samp$SAMPID)
#
samp_subj <- unique(samp[,.(SUBJID,SEX)])
writeLines(sprintf("Total unique subject IDs: %d", length(unique(samp_subj$SUBJID))))
tbl <- table(samp_subj$SEX)
writeLines(sprintf("Subjects (SEX = %s): %d (%.1f%%)", names(tbl), tbl, 100*tbl/nrow(samp_subj)))
#
