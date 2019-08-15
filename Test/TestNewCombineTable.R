library(data.table)
library(foreach)

UpdateCell <- function(current_entry, new_entry) {
  if (length(new_entry) > 0 && !is.na(new_entry)) {
    if (length(new_entry) == 1) {
      # update if no previous entry
      if (length(current_entry) <= 0 || is.na(current_entry)) {
        return(paste0(new_entry, collapse = '.'))
        
        # update if new entry is novel
      } else if (!grepl(new_entry, current_entry, fixed = TRUE, useBytes = TRUE)) {
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
        current_entry <- UpdateCell(current_entry, i)
      }
      return(current_entry)
    }
    
    # do not update if no new entry
  } else if (length(current_entry) > 0 && !is.na(current_entry)) {
    return(current_entry)
    
    # entry remains NA
  } else {
    return(NA)
  }
}

# Given a previously existing entry, add only new information
UpdateEntry <- function(old_entry, entry) {
  new_entry <- old_entry
  if (is.character(entry)) {
    if (is.character(old_entry)) {
      # get list of metadata in new entry
      new_values <- as.list(strsplit(entry, ", ")[[1]])
      
      # check, only add new entries
      for (j in 1:length(new_values)) {
        if (!grepl(new_values[j], new_entry, Encoding("UTF-8"))) {
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
        check.row[1, j] <- UpdateCell(clean.table[row_num, j], table[i, j])
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
          new_entry <- UpdateCell(new_entry, clean.single.trait(j))
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

dict_key <-
  c(
    'Species',
    'Domain',
    'Kingdom',
    'Phylum',
    'Class',
    'Order',
    'Family',
    'Genus',
    'Metabolite Utilization',
    'Enzymes',
    'pH Tolerance',
    'Gram Stain',
    'Energy Source',
    'Growth in Groups',
    'Salt Tolerance',
    'Motility',
    'Oxygen Requirement',
    'Cell Shape',
    'Cell Length',
    'Cell Width',
    'Sporulation',
    'Temperature Range',
    "Type of Hemolysis",
    "Hemolysis Ability",
    "Colony Color",
    "Colony Size",
    "Colony Shape",
    "Incubation Period",
    "Type of Spore",
    "Murein Short Key",
    "Murein Types",
    "Antibiotic Sensitivity",
    "Antibiotic Resistance",
    "Metabolite Production",
    "Pathogenic in Humans",
    "Pathogenic in Animals",
    "Pathogenic in Plants",
    "Flagella",
    'Ecosystem',
    'Pathogenic',
    'Metabolism Type',
    'Mobility',
    'Radioresistance', 
    'Cell Arrangement',
    'Biotic Relationship',
    'Number of membranes'
  )

dict_value <- 
  list(
    list('specie', 'Species', 'species name'),
    list('Domain'),
    list('Kingdom'),
    list('Phylum'),
    list('Class'),
    list('Order'),
    list('Family'),
    list('Genus', 'Genus name'),
    list('metabolite usage', 'Metabolite Utilization', 'Sole carbon substrate use'),
    list('enzymes'),
    list('pH', 'pH range at which growth occurred'),
    list('gram stain', 'Gram stain', 'gram status', 'Gram staining properties'),
    list('energy source', 'Nutrition type', 'Energy source'),
    list('growth in groups', 'Multicellular complex forming ability'),
    list('halophilic', 'halophily', 'NaCl concentration range at which growth occurred (%)'),
    list('motility', 'Motility'),
    list('oxygen requirement', 'Oxygen tolerance', 'oxygen preference', 'Oxygen requirements'),
    list('shape', 'Cell shape', 'cell shape', 'Shape'),
    list('Cell length', 'mean length'),
    list('Cell width', 'mean width'),
    list('sporulation', 'Ability of spore formation', 'spore production', 'Sporulation'),
    list('temperature range', 'Temperature range'),
    list("Type of hemolysis"),
    list("Hemolysis Ability"),
    list("Colony color"),
    list("Colony size"),
    list("Colony shape"),
    list("Incubation period"),
    list('Type of spore'),
    list("Murein short key"),
    list("Murein types"),
    list("Antibiotic Sensitivity"),
    list("Antibiotic Resistance"),
    list("Metabolite Production"),
    list("Pathogenic in humans"),
    list("Pathogenic in animals"),
    list("Pathogenic in plants"),
    list('Flagella', 'flagella', 'Flagellar presence'),
    list('ecosystem', 'Habitat'),
    list('pathogenic', 'Pathogenicity'),
    list('metabolism', 'Metabolism assays', 'Metabolism'),
    list('mobility', 'Mobility'),
    list('radioresistance'),
    list('cell arrangement', 'cell aggregation', 'Cell arrangement'),
    list('biotic relationship', 'Biotic relationship'),
    list('Number of membranes')
  )

names(dict_value) <- dict_key

clean.table <- function(table) {
  # add columns to each data frame so that have same number / names
  for (i in 1:length(table)) {
    names(table)[i] <- gsub('[^a-zA-Z ]+', "", names(table)[i], perl = TRUE)
    
    correct_name <-
      dict_key[which(grepl(
        paste0('\\b',names(table)[i], '\\b'),
        dict_value,
        ignore.case = TRUE
      ))]
    if(length(correct_name) > 0) {
      names(table)[i] <- correct_name
    }
  }
  
  # remove excess
  table <- table[, names(table) %in% dict_key]
  
  # add missing
  for (name in dict_key) {
    if (!(name %in% names(table))) {
      new_column <- data.frame(matrix(NA, nrow = nrow(table), ncol = 1))
      names(new_column) <- name
      table <- cbind(table, new_column)
    }
  }
  
  return(table)
}

protrait <- read.csv('Data/Other/ProTrait_v2019-07-29.csv', stringsAsFactors = FALSE, check.names = FALSE)
ijsem <- read.csv('Data/IJSEM_v2019-08-14.csv', stringsAsFactors = FALSE, check.names = FALSE)
bacmap <- read.csv('Data/Other/BacMap_v2019-08-08.csv',stringsAsFactors = FALSE, check.names = FALSE)
bacdive <- read.csv('Data/Other/BacDive_v2019-07-30.csv', stringsAsFactors = TRUE, check.names = FALSE)

protrait <- clean.table(protrait)
ijsem <- clean.table(ijsem)
bacmap <- clean.table(bacmap)
bacdive <- clean.table(bacdive)

# for (i in 1:nrow(ijsem)) {
#   ijsem$Species[i] <- paste0(ijsem$Genus[i], ' ', ijsem$Species[i])
# }

# rbind all data frames
total_table <- rbindlist(list(protrait, ijsem, bacmap, bacdive), fill = TRUE, use.names = TRUE)

total_table <- clean.repeat(total_table)
message("Merged duplicate species entries")
total_table <- clean.trait(total_table)
message("Removed entry synonyms")
total_table <- clean.nutrition(total_table)
message("Renamed nutrition requirements")
total_table <- apply(total_table, c(1, 2), order.string)
total_table <- as.data.frame(total_table, stringsAsFactors = FALSE, check.names = FALSE)
message("Sorted entries by alphabetical order")

# do clean up methods 

l <- dict_value[which(dict_key == 'a')][[1]][[1]]
l <- dict_key[which(grepl('specie', dict_value, Encoding("UTF-8")))]


for (i in 1:length(total_table)) {
  for (j in 1:nrow(total_table)) {
    if (length(total_table[j, i]) == 0) {
      print(paste0(i, " ", j))
    }
  }
}


# create dictionary from data frame
dictionary <- read.csv('Data/TableDictionary.csv')

