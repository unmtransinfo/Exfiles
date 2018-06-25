#!/usr/bin/env Rscript
#############################################################################
### Copy of Oleg's R, commented.
### GTEx input files all from https://www.gtexportal.org/home/datasets.
### TAU is tissue specificity metric from Yanai et al., 2004,
### https://www.ncbi.nlm.nih.gov/pubmed/15388519.
#############################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)
library(Hmisc)

### Percentile label
perc.label <- function(x, id) {
  tryCatch(
    {
      if(sum(x == 0) == length(x)) {
        return(rep.int(0, length(x)))
      }
      x <- round(x, 3)
      std <- sd(x)
      result <- rep.int(0, length(x))
      q <- quantile(x[which(x >= std)], probs = seq(0,1,0.1))
      if(sum(is.na(q)) == length(q)) {
        return(rep.int(0, length(x)))
      }
      q <- unique(q[!is.na(q)])
      if(length(q) == 1) {
        result[x <= q[1] & x >= std] <- 10
        return(result)
      }
      if(length(q) > 2) {
        for(i in seq(2,length(q) - 1)) {
          if(q[i] == 0 | is.na(q[i])) {
            next
          }
          result[x < q[i] & x >= q[i-1]] <- i - 1
        }
      }
      if(q[length(q)] > 0) {
        result[x <= q[length(q)] & x >= q[length(q) - 1]] <- length(q) - 1
      }
      return(result)
    }, error = function(c) {
      c$message <- paste0(c$message, " (in ", id, ")")
      stop(c)
    }
  )
}

###
tau <- function(x) {
  max.x <- max(x)
  if(max.x == 0) {
    return(0)
  }
  return(sum(sapply(x, function(t) {1-t/max.x}))/(length(x) - 1))
}

###
### READ TISSUE SAMPLES:
### SMTS = tissue name, SMTSD = tissue description
sample.label <- fread("/home/data/GTEx/data/GTEx_v7_Annotations_SampleAttributesDS.txt", header = T, sep = "\t", quote = "", na.strings = "")
### Clean, remove, rename, split cols:
sample.label <- sample.label[, .(SAMPID, SMATSSCR, SMTS, SMTSD)]
sample.label[!is.na(SMTSD),]
# SUBJID is first two hyphen-delimited fields.
sample.label[, c("C1","C2") := tstrsplit(SAMPID, "-", fixed = T, keep = c(1,2))]
sample.label[, SUBJID := sprintf("%s-%s", C1, C2)]
sample.label[, `:=`(C1 = NULL, C2 = NULL)]
#
### READ SUBJECTS:
subject.label <- fread("/home/data/GTEx/data/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", header = T, sep = "\t", quote = "", na.strings = "")
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
sample.label <- sample.label[SMATSSCR < 2 | is.na(SMATSSCR)] # 
sample.label[, `:=`(SUBJID = NULL, SMATSSCR = NULL, SMTS = NULL, SEX = NULL, DTHHRDY = NULL)]
setkey(sample.label, SAMPID)
#
###
### READ GENE TPMs:
### 56202 data rows x 11688 data cols.
### Big: 2.6GB uncompressed.  Need: ~10GB RAM.
### Rows: genes; cols: Name, Description, SampId1, SampId2, ... SampId11688
#rnaseq <- fread("gunzip -c /home/data/GTEx/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz", header = T, skip = 2, sep = "\t", quote = "")
rnaseq <- fread("gunzip -c /home/data/GTEx/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_00001-01000.gct.gz", "\t", skip=2, header=T, sep="\t", quote="") #DEBUG
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
write.table(rnaseq, paste("data/gtex.samples", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = T)
system(paste0("gzip -9v ",paste("data/gtex.samples", Sys.Date(), "tsv", sep = ".")))
###
### Median TPM by gene-tissue:
rnaseq.tau <- rnaseq[, .(MEDIAN_TPM = median(TPM, na.rm = T)), by = .(ENSG, SMTSD)]
###
### TAU = tissue specificity index, by gene (Yanai et al., 2004)
rnaseq.tau <- rnaseq.tau[, .(TAU = tau(perc.label(MEDIAN_TPM, ENSG))), by = ENSG]
#
write.table(rnaseq.tau, paste("data/gtex.tau", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = T)
rm(rnaseq.tau)
system(paste0("gzip -9v ",paste("data/gtex.tau", Sys.Date(), "tsv", sep = ".")))
###
# Assign Low/Medium/High for gene-tissue median rank (percentiles of medians) among tissues.
# Low-Med and Med-High cutoffs are 25 and 75.
# log10(TPM) conventional
setkey(rnaseq, ENSG)
rnaseq.level <- rnaseq[, .(MEDIAN_TPM = median(TPM, na.rm = T)), by = .(ENSG, SMTSD)]
rnaseq.level[, RANK := frank(MEDIAN_TPM)/.N, by = ENSG]
rnaseq.level[MEDIAN_TPM == 0.0, RANK := 0.0]
rnaseq.level[RANK == 0.0, LEVEL := "Not detected"]
rnaseq.level[RANK > 0.0 & RANK < 0.25, LEVEL := "Low"]
rnaseq.level[RANK >= 0.25 & RANK < 0.75, LEVEL := "Medium"]
rnaseq.level[RANK >= 0.75, LEVEL := "High"]
rnaseq.level[, RANK := NULL]
rnaseq.level[, LOG_MEDIAN_TPM := ifelse(MEDIAN_TPM > 0.0, log10(MEDIAN_TPM), NA)]
rnaseq.level[, AGE := "ALL"][, GENDER := "ALL"]
write.table(rnaseq.level, paste("data/gtex.tpm.qualitative", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = T)
rm(rnaseq.level)
#
# Same but by sex.
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
write.table(rnaseq.level.sex, paste("data/gtex.tpm.qualitative", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = F, append = T)
#
# Statify by age.
rnaseq.level.sex.age <- rnaseq[, .(MEDIAN_TPM = median(TPM, na.rm = T)), by = .(ENSG, SMTSD, GENDER, AGE)]
rnaseq.level.sex.age[, RANK := frank(MEDIAN_TPM)/.N, by = ENSG]
rnaseq.level.sex.age[MEDIAN_TPM == 0.0, RANK := 0.0]
rnaseq.level.sex.age[RANK == 0.0, LEVEL := "Not detected"]
rnaseq.level.sex.age[RANK > 0.0 & RANK < 0.25, LEVEL := "Low"]
rnaseq.level.sex.age[RANK >= 0.25 & RANK < 0.75, LEVEL := "Medium"]
rnaseq.level.sex.age[RANK >= 0.75, LEVEL := "High"]
rnaseq.level.sex.age[, RANK := NULL]
rnaseq.level.sex.age[, LOG_MEDIAN_TPM := ifelse(MEDIAN_TPM > 0.0, log10(MEDIAN_TPM), NA)]
setcolorder(rnaseq.level.sex.age, c("ENSG","SMTSD","MEDIAN_TPM","LEVEL","LOG_MEDIAN_TPM","AGE","GENDER"))
write.table(rnaseq.level.sex.age, paste("data/gtex.tpm.qualitative", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = F, append = T)
rm(rnaseq.level.sex.age)
system(paste0("gzip -9v ",paste("data/gtex.tpm.qualitative", Sys.Date(), "tsv", sep = ".")))
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
write.table(sex.diff, paste("data/gtex.sex.diff", Sys.Date(), "tsv", sep = "."), row.names = F, sep = "\t", quote = T, col.names = T)
if(file.exists(paste("data/gtex.sex.diff", Sys.Date(), "tsv.gz", sep = "."))) {
  file.remove(paste("data/gtex.sex.diff", Sys.Date(), "tsv.gz", sep = "."))
}
system(paste0("gzip -9v ",paste("data/gtex.sex.diff", Sys.Date(), "tsv", sep = ".")))
