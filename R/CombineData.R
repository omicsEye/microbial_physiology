library(data.table)
library(foreach)

source('R/Utility.R')

order.string <- function(string) {
  return(paste(sort(unique(trimws(as.vector(strsplit(string, ",")[[1]])))), collapse = ','))
}

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

# create dictionary from data frame
dictionary <- read.csv('Data/TableDictionary.csv', stringsAsFactors = FALSE, check.names = FALSE)
dictionary <- data.frame(lapply(dictionary, trimws), stringsAsFactors = FALSE, check.names = FALSE)

clean.table <- function(table) {
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

combine.data <- function(data, save_file = TRUE) {
  total_table <- list()
  for (i in 1:length(data)) {
    total_table <- rbindlist(list(total_table, clean.table(data[[i]])), fill = TRUE, use.names = TRUE)
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
