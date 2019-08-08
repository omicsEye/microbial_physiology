library(data.table)
library(foreach)

# Given a previously existing entry, add only new information
UpdateEntry <- function(old_entry, entry) {
  new_entry <- old_entry
  if (is.character(entry)) {
    if (is.character(old_entry)) {
      # get list of metadata in new entry
      new_values <- as.list(strsplit(entry, ", ")[[1]])
      
      # check, only add new entries
      for (j in 1:length(new_values)) {
        if (!grepl(new_values[j], new_entry)) {
          new_entry <- paste(new_entry, new_values[j], sep = ", ")
        }
      }
      return(new_entry)
      
      # if no previous entry
    } else {
      return(entry)
    }
    
    # unable to update entry
  } else {
    return(old_entry)
  }
}

order.string <- function(string) {
  return(paste(sort(unique(trimws(as.vector(strsplit(string, ",", useBytes = TRUE)[[1]])))), collapse = ','))
}

clean.repeat <- function(table) {
  clean.table <-
    data.frame(matrix(NA, nrow = 0, ncol = length(table)))
  names(clean.table) <- names(table)
  
  table <- data.frame(lapply(table, as.character), stringsAsFactors=FALSE, check.names = FALSE)
  
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

# Have dictionary where know which columns correspond between tables
# keys are for ProTrait table
dict_key <- c(
  "specie",
  "metabolite usage",
  "enzymes",
  "pH",
  "gram stain",
  "energy source",
  "growth in groups",
  "halophilic",
  "motility",
  "oxygen requirement",
  "shape",
  "sporulation",
  "temperature range"
)

# values are for BacDive table
dict_value <- c(
  "Species",
  "Metabolite Utilization",
  "enzymes",
  "pH",
  "Gram stain",
  "Nutrition type",
  "Multicellular complex forming ability",
  "halophily",
  "Motility",
  "Oxygen tolerance",
  "Cell shape",
  "Ability of spore formation",
  "Temperature range"
)

# Create dictionary
names(dict_value) <- dict_key

# Trait labels from ProTraits
protrait_names <-
  c(
    "specie",
    "metabolite usage",
    "enzymes",
    "ecosystem",
    "pH",
    "biotic relationship",
    "cell arrangement",
    "energy source",
    "flagella",
    "gram stain",
    "growth in groups",
    "halophilic",
    "pathogenic",
    "metabolism",
    "mobility",
    "motility",
    "oxygen requirement",
    "radioresistance",
    "shape",
    "sporulation",
    "temperature range"
  )

# Trait lavels from BacDive
bacdive_names <- c(
  "Specie",
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Gram stain",
  "Cell length",
  "Cell width",
  "Cell shape",
  "Motility",
  "Type of hemolysis",
  "Hemolysis Ability",
  "Colony color",
  "Colony size",
  "Colony shape",
  "Incubation period",
  "Ability of spore formation",
  "Type of spore",
  "Multicellular complex forming ability",
  "Name of produced compound",
  "Murein short key",
  "Murein types",
  "Oxygen tolerance",
  "Nutrition type",
  "Antibiotic Sensitivity",
  "Antibiotic Resistance",
  "halophily",
  "Metabolite Utilization",
  "Metabolite Production",
  "enzymes",
  "Temperature range",
  "pH",
  "Pathogenic in humans",
  "Pathogenic in animals",
  "Pathogenic in plants"
)

# For testing purposes
# protrait <- read.delim(file.choose(), fill=NA, stringsAsFactors = FALSE, sep=',')
# names(protrait) <- protrait_names
# bacdive <- read.delim(file.choose(), fill=NA, stringsAsFactors = FALSE, sep=',')
# names(bacdive) <- bacdive_names

# combine based on common species
# combine protrait table to bacdive table
CombineData <- function(protrait, bacdive, save_file = TRUE) {
  ctable <- bacdive  # new table that will contain combined table
  
  # add columns to table that does not have already (in protraits, not bacdive)
  for (name in protrait_names) {
    if (!(name %in% dict_key)) {
      new_column <- data.frame(matrix(NA, nrow = nrow(ctable), ncol = 1))
      names(new_column) <- name
      ctable <- cbind(ctable, new_column)
    }
  }
  
  # add information from protraits to new table
  num_pro_rows <- nrow(protrait)
  num_pro_cols <- ncol(protrait)
  num_new_cols <- ncol(ctable)
  new_row_names <- names(ctable) # column names for new table
  
  row_count <- 1
  while (row_count <= num_pro_rows) {
    # combine by microbe
    current_pro_row <- protrait[row_count,]
    new_row <- data.frame(matrix(NA, nrow = 1, ncol = num_new_cols))
    names(new_row) <- new_row_names
    
    col_count <- 1
    # check if specie already exists
    match <- which(ctable[, 1] %in% new_row[1])
    if (length(match) > 0) {
      # update first instance
      update_row <- ctable[match[1], ]
      
      while (col_count <= num_pro_cols) {
        cell <- current_pro_row[col_count]
        column_name <- names(current_pro_row)[col_count]
        
        # translate label if necessary
        if (column_name %in% dict_key) {
          column_name <- dict_value[column_name]
        }
        
        # add trait to correct cell
        update_row[column_name] <-
          UpdateEntry(update_row[column_name], cell)
        
        col_count <- col_count + 1
      }
      
      ctable[match[1], ] <- update_row
      
    } else {
      while (col_count <= num_pro_cols) {
        cell <- current_pro_row[col_count]
        column_name <- names(current_pro_row)[col_count]
        
        # translate label if necessary
        if (column_name %in% dict_key) {
          column_name <- dict_value[column_name]
        }
        
        # add trait to correct cell
        new_row[column_name] <- cell
        
        col_count <- col_count + 1
      }
      
      # update table
      ctable <- rbind(ctable, new_row)
    }
    # move on to next microbe
    row_count <- row_count + 1
  }
  
  message("Completed table merge")
  
  ctable <- clean.repeat(ctable)
  message("Merged duplicate species entries")
  ctable <- clean.trait(ctable)
  message("Removed entry synonyms")
  ctable <- clean.nutrition(ctable)
  message("Renamed nutrition requirements")
  ctable <- apply(ctable, c(1, 2), order.string)
  ctable <- as.data.frame(ctable, stringsAsFactors = FALSE, check.names = FALSE)
  message("Sorted entries by alphabetical order")
  
  if (save_file) {
    write.csv(ctable,
              paste0("mxp_microbiome_v", Sys.Date(), ".csv"),
              row.names = FALSE)
  }
  
  return(ctable)
}
