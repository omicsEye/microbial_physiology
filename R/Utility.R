library(foreach)
library(RAM)
library(recommenderlab)
library(readxl)

#' returns data frame of excel file
#'
#' @param path the path to the excel file
#'
#' @return a dataframe representation
#' @export
read.excel <- function(path) {
  return(as.data.frame(read_excel(path)))
}

#' Remove all entries that do not appear enough
#'
#' @param table table to clean
#' @param percentage threshold to eliminate entries
#'
#' @return cleaned table
#' @export
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

#' load and prepare abundance data from .csv
#' @details recommeneded for metabolite abundance tables
#' @param path path to abundance table
#' @param column
#' @return a numerical matrix of abundance data
#' @export
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

#' load and prepare meta data from .csv
#' @param path path to the meta data
#' @param tax_column number of column containing naming infomration
#' @return a data frame
#' @export
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

#' load and prepare ihmp metabolite intensity value table from .csv
#' @details specialized for opening ihmp abundance table files
#'
#' @param file path to file containg ihmp data
#'
#' @return a data frame
#' @export
load.ihmp <- function(file) {
  ihmp <- read.csv(
    file,
    header = TRUE,
    fill = TRUE,
    comment.char = "" ,
    check.names = TRUE,
    stringsAsFactors = FALSE
  )
  ihmp_col_names <- ihmp[4,-c(1:7)]

  # filter out for only hmdb id
  ihmp  <- ihmp[(ihmp$X.5) != "",]

  # get column_names
  ihmp_row_names <- ihmp[-1,6]

  ihmp <- ihmp[-1,-c(1:7)]

  ihmp <- apply(as.matrix(ihmp), 2, as.numeric)

  # clean up ihmp data
  colnames(ihmp) <- ihmp_col_names
  rownames(ihmp) <- ihmp_row_names
  ihmp[is.na(ihmp)] <- 0
  return(ihmp)
}

#' normalize ihmp metabolite intensity values using z-score normalization
#'
#' @param ihmp ihmp matrix
#'
#' @return normalized matrix
#'
normalize.ihmp <- function(ihmp) {
  ihmp[is.na(ihmp)] <- 0
  ihmp <- as(ihmp, "realRatingMatrix")
  ihmp <- normalize(ihmp, method="Z-score", row=TRUE)
  ihmp <- as.matrix(ihmp@data)
  ihmp[is.na(ihmp)] <- 0
  return(ihmp)
}

#' given a cell entry, update it based on a new entry
#'
#' @param current_entry information containined in current cell
#' @param new_entry information that want to add to cell
#'
#' @return cell with updated information
update.cell <- function(current_entry, new_entry) {
  if (length(new_entry) > 0 && !is.na(new_entry)) {
    if (length(new_entry) == 1) {
      # update if no previous entry
      if (length(current_entry) <= 0 || is.na(current_entry)) {
        return(paste0(new_entry, collapse = '.'))

        # update if new entry is novel
      } else if (!grepl(new_entry, current_entry, fixed = TRUE)) {
        return(paste(
          current_entry,
          new_entry,
          sep = ', ',
          collapse = ", "
        ))

        # do not update if new entry not novel
      } else {
        return(current_entry)
      }

      # split if mutliple items in entry
    } else if (length(new_entry) >= 1) {
      for (i in new_entry) {
        current_entry <- update.cell(current_entry, i)
      }
      return(current_entry)
    }

    # do not update if no new entry
  } else if (!is.na(current_entry)) {
    return(current_entry)

    # entry remains NA
  } else {
    return(NA)
  }
}
