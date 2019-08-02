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
