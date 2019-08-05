library(foreach)

# Remove all entries that do not appear enough
clear.small.entry <- function(table, percentage) {
  for(i in 1:ncol(table)) {
    terms <- unique(table[, i])
    for(j in 1:length(terms)) {
      if(length(which(table[, i] == terms[j])) / length(table[, i]) < percentage) {
        table[table == terms[j]] <- NA
      }
    }
  }
  table[table == ''] <- NA
  return(table)
}

# load and prepare abundance data from .csv
# returns a numerical matrix
load.abundance.data <- function(path, column = 1) {
  if (is.na(path)) {
    path <- file.choose()
  }
  abundance_table <-
    read.delim(
      path,
      sep = ",",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  abundance_table <-
    abundance_table[!is.na(abundance_table[, column]),]
  abundance_table <-
    abundance_table[(abundance_table[, column]) != '',]
  abundance_table <-
    abundance_table[!duplicated(abundance_table[, column]),]
  species <- abundance_table[, column]
  abundance_table <- apply(abundance_table, 2, as.numeric)
  row.names(abundance_table) <- species
  return(abundance_table)
}

# load and prepare meta data from .csv
# returns a data frame
load.meta.data <- function(path, tax_column = 1) {
  if (is.na(path)) {
    path <- file.choose()
  }
  data <-
    read.delim(
      path,
      fill = NA,
      stringsAsFactors = FALSE,
      sep = ',',
      check.names = FALSE
    )
  data <- data[!is.na(data[, tax_column]),]
  data <- data[(data[, tax_column]) != '',]
  data <- data[!duplicated(data[, tax_column]),]
  row.names(data) <- data[, tax_column]
  data <- data[,-c(tax_column)]
  return(data)
}
