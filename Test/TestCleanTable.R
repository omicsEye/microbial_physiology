library(data.table)
library(foreach)

order.string <- function(string) {
  return(paste(sort(unique(trimws(as.vector(strsplit(string, ",")[[1]])))), collapse = ','))
}

clean.repeat <- function(table) {
  clean.table <-
    data.frame(matrix(NA, nrow = 0, ncol = length(table)))
  names(clean.table) <- names(table)
  
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE)
  
  for(i in 1:nrow(table)) {
    
    if(table$Species[i] %in% clean.table$Species) {
      row_num <- which(clean.table$Species == table$Species[i])
      check.row <- clean.table[row_num, , drop = FALSE]
      
      for(j in 1:ncol(check.row)) {
        check.row[1, j] <- UpdateCell(clean.table[row_num, j], table[i, j])
      }
      clean.table <- rbind(clean.table[-c(row_num), ], check.row, stringsAsFactors = FALSE)
      
    } else {
      clean.table <- rbind(clean.table, table[i, ], stringsAsFactors = FALSE)
    }
  }
  return(clean.table)
}

clean.nutrition <- function(table) {
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)
  
  foreach(i = 1:nrow(table)) %do% {
    new_entry <- ""
    current_entry <- table$`Nutrition type`[i]
    
    if(!is.na(current_entry)) {
      # check each nutrition type
      if(grepl("auto", current_entry))
        new_entry <- paste0(new_entry, "autotroph, ")
      if(grepl("litho", current_entry))
        new_entry <- paste0(new_entry, "lithotroph, ")
      if(grepl("chemo", current_entry))
        new_entry <- paste0(new_entry, "chemotroph, ")
      if(grepl("hetero", current_entry))
        new_entry <- paste0(new_entry, "heterotroph, ")
      if(grepl("organo", current_entry))
        new_entry <- paste0(new_entry, "organotroph, ")
      if(grepl("methylo", current_entry))
        new_entry <- paste0(new_entry, "methylotroph, ")
      if(grepl("photo", current_entry))
        new_entry <- paste0(new_entry, "phototroph, ")
      
      table$`Nutrition type`[i] <- order.string(substring(new_entry, 1, nchar(new_entry) - 2))
    }
  }
  return(table)
}

clean.trait <- function(table) {
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)
  
  foreach(i = 1:nrow(table)) %do% {
    trait <- table$`Oxygen tolerance`[i]
    
    if(!is.na(trait)) {
      individual <- trimws((unlist(strsplit(trait, ","))))
      
      if(length(individual) == 1) {
        table$`Oxygen tolerance`[i] <- clean.single.trait(trait)
        
      } else {
        new_entry <- ''
        for(j in individual) {
          new_entry <- UpdateCell(new_entry, clean.single.trait(j))
        }
        table$`Oxygen tolerance`[i] <- order.string(substring(new_entry, 3))
      }
    }
  }
  return(table)
}

clean.single.trait <- function(trait) {
  trait_dictionary <-
    list(
      list("anaerobe", "strictanaero", "anaero", "obligate anaerobe"),
      list("aerobe", "strictaero", "aero", "obligate aerobe"),
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
