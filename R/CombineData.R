library(data.table)
library(foreach)

#' Order string
#' Order all words in a string, which are separated by commas, in order
#'
#' @param string the string to order
#'
#' @return the string in alphabetical order
order.string <- function(string) {
  return(paste(sort(unique(trimws(as.vector(strsplit(string, ",")[[1]])))), collapse = ','))
}

#' Remove repeat entries
#' remove all entries that have the same species name and merge all unique data
#'
#' @param table the table to clean
#'
#' @return the cleaned table
clean.repeat <- function(table) {
  clean.table <-
    data.frame(matrix(NA, nrow = 0, ncol = length(table)))
  names(clean.table) <- names(table)

  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)

  for(i in 1:nrow(table)) {

    if(!is.na(table$Species[i]) && table$Species[i] %in% clean.table$Species) {
      row_num <- which(clean.table$Species == table$Species[i])
      check.row <- clean.table[row_num, , drop = FALSE]

      for(j in 1:ncol(check.row)) {
        check.row[1, j] <- update.cell(clean.table[row_num, j], table[i, j])
      }
      clean.table <- rbind(clean.table[-c(row_num), ], check.row, stringsAsFactors = FALSE)

    } else {
      clean.table <- rbind(clean.table, table[i, ], stringsAsFactors = FALSE)
    }
  }
  return(clean.table)
}

#' split all nutrition data types into its individual components
#'
#' @param table the table containing nutrition requirement information
#' @param trait the name of the column containing nutrition requirement information
#'
#' @return the cleaned table
clean.nutrition <- function(table, trait = 'Energy Source') {
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)

  foreach(i = 1:nrow(table)) %do% {
    new_entry <- ""
    current_entry <- table[i, trait]

    if(!is.na(current_entry)) {
      # check each nutrition type
      if(grepl("auto", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "autotroph, ")
      if(grepl("litho", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "lithotroph, ")
      if(grepl("chemo", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "chemotroph, ")
      if(grepl("hetero", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "heterotroph, ")
      if(grepl("organo", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "organotroph, ")
      if(grepl("methylo", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "methylotroph, ")
      if(grepl("photo", current_entry, Encoding("UTF-8")))
        new_entry <- paste0(new_entry, "phototroph, ")

      table[i, trait] <- order.string(substring(new_entry, 1, nchar(new_entry) - 2))
    }
  }
  return(table)
}

#' clean trait entries that have synonyms
#'
#' @param table the table to be cleaned
#' @param trait_name the name of the column that needs to be cleaned
#'
#' @return the cleaned table
clean.trait <- function(table, trait_name = 'Oxygen Requirement') {
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)

  foreach(i = 1:nrow(table)) %do% {
    trait <- table[i, trait_name]

    if(!is.na(trait)) {
      individual <- trimws((unlist(strsplit(trait, ","))))

      if(length(individual) == 1) {
        table[i, trait_name] <- clean.single.trait(trait)

      } else {
        new_entry <- ''
        for(j in individual) {
          new_entry <- update.cell(new_entry, clean.single.trait(j))
        }
        table[i, trait_name] <- order.string(substring(new_entry, 3))
      }
    }
  }
  return(table)
}

#' clean the trait column (currently only for oxygen requirements)
#' @param trait the synoymn
#'
#' @return the predetermined name based on the synonym given
clean.single.trait <- function(trait) {
  trait_dictionary <-
    list(
      list("anaerobe", "strictanaero", "anaero", "obligate anaerobe", "anaerobic"),
      list("aerobe", "strictaero", "aero", "obligate aerobe", "aerobic"),
      list("facultative aerobe", "facultative anaerobe", "facultative"),
      list("microaerophile", "microaerophilic")
    )
  names(trait_dictionary) <- c("anaerobe", "aerobe", "facultative", "microaerophile")

  for(j in 1:length(trait_dictionary)) {
    if(any(trait == trait_dictionary[[j]])) {
      return(names(trait_dictionary)[j])
    }
  }
  return(trait)
}

#' create dictionary from data frame
#'
#' @param path the path to a csv containg trait names and synonyms
#'
#' @return the data frame representation of the csv
read.trait.dictionary <- function(path) {
  dictionary <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  dictionary <- data.frame(lapply(dictionary, trimws), stringsAsFactors = FALSE, check.names = FALSE)
  return(dictionary)
}

#' prepare table so that can be combined with other tables
#'
#' @param table the table to be prepared
#' @param path path to a table contianing trait names and synoynms
#'
#' @return the prepared table
clean.table <- function(table, path) {
  dictionary <- read.trait.dictionary(path)

  # add columns to each data frame so that have same number / names
  for (i in 1:length(table)) {
    names(table)[i] <- trimws(gsub('[^a-zA-Z ]+', "", names(table)[i], perl = TRUE))

    correct_name <-
      names(dictionary)[which(names(table)[i] == dictionary, arr.ind = TRUE)[2]]
    print(which(names(table)[i] == dictionary, arr.ind = TRUE))
    if(length(correct_name) > 0) {
      names(table)[i] <- correct_name
    }
  }

  # remove excess
  table <- table[, names(table) %in% names(dictionary)]

  # add missing
  for (name in names(dictionary)) {
    if (!(name %in% names(table))) {
      new_column <- data.frame(matrix(NA, nrow = nrow(table), ncol = 1))
      names(new_column) <- name
      table <- cbind(table, new_column)
    }
  }

  return(table)
}

#' Combine data
#' Combine a list of tables together
#' @details the csv that must be supplied needs to contain the original column names and its synonym or it will be removed
#' @param data a list of data frames to combine
#' @param path to a table contianing trait names and synoynms
#'
#' @return a single data frame containing the combined and clean up information
#' @export
combine.data <- function(data, path = 'Data/TableDictionary.csv', save_file = TRUE) {
  total_table <- list()
  for (i in 1:length(data)) {
    total_table <- rbindlist(list(total_table, clean.table(data[[i]], path)), fill = TRUE, use.names = TRUE)
  }

  total_table <- clean.repeat(total_table)
  message("Merged duplicate species entries")
  total_table <- clean.trait(total_table)
  message("Removed entry synonyms")
  total_table <- clean.nutrition(total_table)
  message("Renamed nutrition requirements")
  total_table <- apply(total_table, c(1, 2), order.string)
  total_table <- as.data.frame(total_table, stringsAsFactors = FALSE, check.names = FALSE)
  message("Sorted entries by alphabetical order")

  if (save_file) {
    write.csv(total_table,
              paste0("mxp_microbiome_v", Sys.Date(), ".csv"),
              row.names = FALSE)
  }

  return(total_table)
}
