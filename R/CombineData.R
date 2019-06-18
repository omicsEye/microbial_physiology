# Given a previously existing entry, add only new information 
UpdateEntry <- function(old_entry, entry) {
  new_entry <- old_entry
  if(is.character(entry)) {
    if(is.character(old_entry)) {
      new_values <- as.list(strsplit(entry, ", ")[[1]])
      
      for(j in 1:length(new_values)) {
        if(!grepl(new_values[j], new_entry)) {
          new_entry <- paste(new_entry, new_values[j], sep=", ")
        }
      }
      return(new_entry)
    } else {
      return(entry)
    }
  } else {
    return(old_entry)
  }
}

# Have dictionary where know which columns correspond between tables
# keys are for ProTrait table
dict_key <- c("specie", "metabolite usage", "enzymes", "pH", 
              "gram stain", "energy source",
              "growth in groups", "halophilic", 
              "motility", 
              "oxygen requirement", "shape", "sporulation", "temperature range")
# values are for BacDive table
dict_value <- c("Species", "met_util", "enzymes", "pH",
                "Gram stain", "Nutrition type",
                "Multicellular complex forming ability", "halophily",
                "Motility", "Oxygen tolerance", "Cell shape", "Ability of spore formation",
                "Temperature range")

names(dict_value) <- dict_key

protrait_names <- c("specie", "metabolite usage", "enzymes", "ecosystem", "pH", 
                    "biotic relationship", "cell arrangement",
                    "energy source", "flagella", "gram stain", "growth in groups", "halophilic", 
                    "pathogenic", "metabolism", "mobility", "motility", 
                    "oxygen requirement", "radioresistance", "shape", "sporulation", "temperature range")

bacdive_names <- c("Specie", "Domain", "Phylum", "Class", "Order",
                   "Family", "Genus",
                   "Gram stain", "Cell length", "Cell width",
                   "Cell shape", "Motility",
                   "Type of hemolysis", "Hemolysis Ability", "Colony color",
                   "Colony size", "Colony shape",
                   "Incubation period",
                   "Ability of spore formation",
                   "Type of spore",
                   "Multicellular complex forming ability",
                   "Name of produced compound",
                   "Murein short key", "Murein types",
                   "Oxygen tolerance",
                   "Nutrition type",
                   "met_antibiotica", "halophily",
                   "met_util", "met_production", "enzymes",
                   "Temperature range", "pH",
                   "Pathogenic in humans", "Pathogenic in animals", 
                   "Pathogenic in plants")

# protrait <- read.delim(file.choose(), fill=NA, stringsAsFactors = FALSE, sep=',')
# names(protrait) <- protrait_names
# bacdive <- read.delim(file.choose(), fill=NA, stringsAsFactors = FALSE, sep=',')
# names(bacdive) <- bacdive_names

# combine based on common species 
# combine protrait table to bacdive table
CombineData <- function(protrait, bacdive, save_file=TRUE) {
  ctable <- bacdive
  
  # add columns to table that does not have already (in protraits, not bacdive)
  for(name in protrait_names) {
    if(!(name %in% dict_key)) {
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
  while(row_count <= num_pro_rows) {
    current_pro_row <- protrait[row_count, ]
    new_row <- data.frame(matrix(NA, nrow = 1, ncol = num_new_cols))
    names(new_row) <- new_row_names
    
    col_count <- 1
    # check if specie already exists 
    if(new_row[1] %in% ctable[ ,1]) {
      
    } else {
      while(col_count <= num_pro_cols) {
        cell <- current_pro_row[col_count]
        column_name <- names(current_pro_row)[col_count]
        if(column_name %in% dict_key) {
          column_name <- dict_value[column_name]
        }
        new_row[column_name] <- cell
        
        col_count <- col_count + 1
      }
      ctable <- rbind(ctable, new_row)
    }
    row_count <- row_count + 1
  }
  if(save_file) {
    write.csv(ctable, "CombineTable.csv", row.names = FALSE)
  }
  message("Completed table merge")
  return(ctable)
}
