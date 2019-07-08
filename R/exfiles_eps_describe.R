library(readr)

eps <- read_delim("data/exfiles_eps.tsv", "\t", col_types=cols(.default=col_number(), ENSG=col_character(), SEX=col_character()))

for (coltag in names(eps)) {
  writeLines(sprintf("N_NA = %4d (%4.1f%%)\t%s", sum(is.na(eps[[coltag]])), 100*sum(is.na(eps[[coltag]]))/nrow(eps), coltag))
}

