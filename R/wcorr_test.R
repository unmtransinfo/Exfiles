#!/usr/bin/env Rscript
#############################################################################
#############################################################################
library(data.table)
library(wCorr)
#
###
wPearson <- function(A, B) { #Vector version (slow)
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
# Corresponding rows of each matrix compared.
wPearson_mx <- function(A, B) { #Matrix version
  mapply(weightedCorr, as.list(as.data.frame(t(A))), as.list(as.data.frame(t(B))), 
         weights=as.list(as.data.frame(t((A+B)/2))),
         MoreArgs=list(method="Pearson"))
}


epA <- rnorm(5, mean=10, sd=2)
epB <- rnorm(5, mean=10, sd=2)

print(wPearson(epA, epB))

#

N_IDCOLS <- 2

epsA <- data.table(ENSG = sprintf("FOOBAR%02d", 1:12), SEX = c(rep("F", 6), rep("M", 6)))
for (j in 3:11) {
  epsA[[sprintf("TISSUE%02d", j)]] <- rnorm(12, mean=0)
}

epsB <- data.table(ENSG = sprintf("FOOBAR%02d", 1:12), SEX = c(rep("F", 6), rep("M", 6)))
for (j in 3:11) {
  epsB[[sprintf("TISSUE%02d", j)]] <- rnorm(12, mean=0)
}

epsA_mx <- as.matrix(epsA[, (N_IDCOLS+1):ncol(epsA)])
epsB_mx <- as.matrix(epsB[, (N_IDCOLS+1):ncol(epsB)])

print(wPearson_mx(epsA_mx, epsB_mx))

for (i in 1:nrow(epsA)) {
  message(sprintf("%d. %f", i, wPearson(epsA_mx[i,], epsB_mx[i,])))
}

