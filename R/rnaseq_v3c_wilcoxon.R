#!/usr/bin/env Rscript
#############################################################################
### GTEx input files all from https://www.gtexportal.org/home/datasets.
#############################################################################
library(readr)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
#library(Hmisc)

fdir <- "/home/data/GTEx/data"
###
### READ TISSUE SAMPLES:
### SMTS = tissue name, SMTSD = tissue description
fname <- "GTEx_v7_Annotations_SampleAttributesDS.txt"
writeLines(sprintf("File: %s", fname))
sample.label <- fread(sprintf("%s/%s", fdir, fname), header = T, sep = "\t", quote = "", na.strings = "")
### Clean, remove, rename, split cols:
sample.label <- sample.label[, .(SAMPID, SMATSSCR, SMTS, SMTSD)]
sample.label[!is.na(SMTSD),]
# SUBJID is first two hyphen-delimited fields.
sample.label[, c("C1","C2") := tstrsplit(SAMPID, "-", fixed = T, keep = c(1,2))]
sample.label[, SUBJID := sprintf("%s-%s", C1, C2)]
sample.label[, `:=`(C1 = NULL, C2 = NULL)]
#
### READ SUBJECTS:
fname <- "GTEx_v7_Annotations_SubjectPhenotypesDS.txt"
writeLines(sprintf("File: %s", fname))
subject.label <- fread(sprintf("%s/%s", fdir, fname), header = T, sep = "\t", quote = "", na.strings = "")
### Keep only subjects healthy at death: 
### (DTHHRDY = 4-point Hardy Scale Death Classification.)
subject.label <- subject.label[DTHHRDY == 1 | DTHHRDY == 2]
#
### MERGE SAMPLES + SUBJECTS
sample.label <- merge(sample.label, subject.label, by.x = "SUBJID", by.y = "SUBJID")
###  Clean, remove, rename cols:
sample.label[, GENDER := factor(SEX, levels = c(1,2), labels = c("male", "female"))]
sample.label[is.na(SMTS) & SMTSD %like% "Skin - ", SMTS := "Skin"]
### Remove samples with high degree of autolysis:
sample.label <- sample.label[SMATSSCR < 2 | is.na(SMATSSCR)]
sample.label[, `:=`(SUBJID = NULL, SMATSSCR = NULL, SMTS = NULL, SEX = NULL, DTHHRDY = NULL)]
setkey(sample.label, SAMPID)
#
###
### READ GENE TPMs:
### 56202 data rows x 11688 data cols.
### Big: 2.6GB uncompressed.  Need: ~10GB RAM.
### Rows: genes; cols: Name, Description, SampId1, SampId2, ... SampId11688
#fname <- "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_00001-01000.gct.gz"
fname <- "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
writeLines(sprintf("File: %s", fname))
rnaseq <- fread(cmd=sprintf("gunzip -c %s/%s", fdir, fname), skip=2, header=T, sep="\t", quote="")
### Clean, remove, rename cols:
rnaseq[, Description:=NULL]
setnames(rnaseq, old = c("Name"), new = c("ENSG"))
### MELT: One row per ENSG-SAMPID-TPM triplet:
rnaseq <- melt(rnaseq, id.vars = "ENSG", variable.name = "SAMPID", value.name = "TPM", na.rm = T, value.factor = F, verbose = T)
setkey(rnaseq, SAMPID)
#
# Remove genes in pseudoautosomal regions (PAR) of chromosome Y ("ENSGR").
rnaseq <- rnaseq[!ENSG %like% "ENSGR"]
### MERGE with samples.
rnaseq <- merge(rnaseq, sample.label, by.x = "SAMPID", by.y = "SAMPID")
setkey(rnaseq, ENSG)
###
# Assign Low/Medium/High for gene-tissue median rank (percentiles of medians) among tissues, by SEX.
# Low-Med and Med-High cutoffs are 25 and 75.
# log10(TPM) conventional
###
rnaseq.level.sex <- rnaseq[, .(MEDIAN_TPM = median(TPM, na.rm = T)), by = .(ENSG, SMTSD, GENDER)]
rnaseq.level.sex[, RANK := frank(MEDIAN_TPM)/.N, by = ENSG]
rnaseq.level.sex[MEDIAN_TPM == 0.0, RANK := 0.0]
rnaseq.level.sex[RANK == 0.0, LEVEL := "Not detected"][RANK > 0.0 & RANK < 0.25, LEVEL := "Low"]
rnaseq.level.sex[RANK >= 0.25 & RANK < 0.75, LEVEL := "Medium"]
rnaseq.level.sex[RANK >= 0.75, LEVEL := "High"]
rnaseq.level.sex[, RANK := NULL]
rnaseq.level.sex[, LOG_MEDIAN_TPM := ifelse(MEDIAN_TPM > 0.0, log10(MEDIAN_TPM), NA)]
rnaseq.level.sex[, AGE := "ALL"]
setcolorder(rnaseq.level.sex, c("ENSG","SMTSD","MEDIAN_TPM","LEVEL","LOG_MEDIAN_TPM","AGE","GENDER"))
#
# Remove data for gene-tissue pairs not present in both sexes.
rnaseq[, sex.count := uniqueN(GENDER), by = .(ENSG, SMTSD)]
rnaseq <- rnaseq[sex.count == 2]
rnaseq[, sex.count := NULL]
#
# Remove data for gene-tissue pairs with all zero expression.
rnaseq[, max.tpm.0 := ifelse(max(TPM) == 0, T, F), by = .(ENSG, SMTSD)]
rnaseq <- rnaseq[max.tpm.0 == F]
rnaseq[, max.tpm.0 := NULL]
#
###
# For each gene-tissue, compute sex difference via Wilcox test.
sex.diff <- rnaseq[, .(p.value = wilcox.test(TPM ~ GENDER, data = .SD)$p.value, statistic = wilcox.test(TPM ~ GENDER, data = .SD)$statistic), by = .(ENSG, SMTSD)]
#
#
gtex <- rnaseq.level.sex
gtex <- gtex[, .(ENSG, SMTSD, MEDIAN_TPM, GENDER)]
#
# Combine rows into one row per gene+tissue, cols for M and F TPM.
gtex <- dcast(gtex, ENSG + SMTSD ~ GENDER, value.var = "MEDIAN_TPM", fun.aggregate = mean, fill = NA, drop = T)
gtex <- gtex[!is.na(female) & !is.na(male)]
#
# Log fold-change is log of ratio.
gtex[, log2.fold.change := log2(max(male,female)/min(male,female)), by = .(ENSG, SMTSD)]
#
gtex[, key := sprintf("%s.%s", ENSG, SMTSD)]
sex.diff[, key := sprintf("%s.%s", ENSG, SMTSD)]
#
sex.diff <- merge(sex.diff, gtex[, .(female, male, log2.fold.change, key)], by.x = "key", by.y = "key")
setnames(sex.diff, c("female","male"), c("MEDIAN_TPM_F", "MEDIAN_TPM_M"))
sex.diff[, key := NULL]
#
write_delim(sex.diff, "data/gtex.sex.diff.tsv", "\t")
#
