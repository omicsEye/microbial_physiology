# adding an entry to a cell so that do not overwrite previous information
AddToCell <- function(cell, new_entry) {
  # if no previous entry
  if (cell == "" | is.na(cell)) {
    return(new_entry)
    # data set known to not have duplicates
  } else {
    return(paste(cell, new_entry, sep = ", "))
  }
}

# get trait from column name given cell
GetTrait <- function(col_names, current_column) {
  trait <- col_names[current_column]
  trait <-
    gsub('[[:punct:]]', ' ', substring(trait, regexpr("\\.", trait)[[1]][1] + 1))
  return(trait)
}

# Organize ProTrait table into clean format
CleanProTrait <- function(save_file = TRUE) {
  # known traits that are going to be extracted
  column_names <-
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
  
  # check if local file exists, download if does not
  if (!file.exists('ProTraits_binaryIntegratedPr0.90.txt')) {
    download.file(
      'http://protraits.irb.hr/data/ProTraits_binaryIntegratedPr0.90.txt',
      'ProTraits_binaryIntegratedPr0.90.txt'
    )
  }
  raw_table <-
    read.delim(
      'ProTraits_binaryIntegratedPr0.90.txt',
      fill = NA,
      stringsAsFactors = FALSE
    )
  # only want data that can be parsed
  raw_table <- raw_table[, 1:313]
  
  # dimensions of raw table
  num_rows <- nrow(raw_table)
  num_cols <- ncol(raw_table)
  
  # new table that will contain all the extracted data
  refined_table <-
    data.frame(matrix(nrow = num_rows, ncol = length(column_names)))
  names(refined_table) <- column_names
  
  current_line <- 1
  
  message("Began ProTrait extraction")
  while (current_line <= num_rows) {
    # extract data one microbe at a time
    current_row <- raw_table[current_line, ]
    
    # prompt user about progress
    if (current_line %% 500 == 0) {
      message(paste("Cleaned", current_line, "/", num_rows, "rows ..."))
    }
    
    # get specie name
    refined_table[current_line, "specie"] <- current_row[1]
    
    # relavent information begins in third column
    current_column <- 3
    
    # extract data one column at a time
    while (current_column <= num_cols) {
      current_cell <- current_row[current_column]
      current_name <- names(current_cell)
      
      # get gram stain - unique because 0 represents gram negative rather than lack of phenotype
      if (grepl("gram_stain", current_name) & !is.na(current_cell)) {
        if (current_cell == 1) {
          refined_table[current_line, "gram stain"] <-
            "positive"
          
        } else if (current_cell == 0) {
          refined_table[current_line, "gram stain"] <-
            "negative"
        }
        # Only want to record traits that microbe has
      } else if (!is.na(current_cell) & current_cell == 1) {
        # get all metabolites and enzymes from column 3 to 111
        if (current_column >= 3 & current_column <= 111) {
          # Find the enzymes
          if (grepl("ase", current_name) |
              grepl("trypsin", current_name)) {
            enzyme_cell <- refined_table[current_line, "enzymes"]
            refined_table[current_line, "enzymes"] <-
              AddToCell(enzyme_cell, GetTrait(names(raw_table), current_column))
            
            # or its a metabolite
          } else {
            meta_cell <- refined_table[current_line, "metabolite usage"]
            refined_table[current_line, "metabolite usage"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get remaining metadata
        } else {
          # get ecosystem
          if (grepl("ecosystem", current_name) |
              grepl("habitat", current_name)) {
            meta_cell <- refined_table[current_line, "ecosystem"]
            refined_table[current_line, "ecosystem"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get pH type
          if (grepl("phenotype", current_name)) {
            meta_cell <- refined_table[current_line, "pH"]
            refined_table[current_line, "pH"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get biotic relationship
          if (grepl("bioticrelationship", current_name)) {
            meta_cell <- refined_table[current_line, "biotic relationship"]
            refined_table[current_line, "biotic relationship"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get cell arrangement
          if (grepl("cellarrangement", current_name)) {
            meta_cell <- refined_table[current_line, "cell arrangement"]
            refined_table[current_line, "cell arrangement"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get energy source
          if (grepl("energysource", current_name)) {
            meta_cell <- refined_table[current_line, "energy source"]
            refined_table[current_line, "energy source"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get flagella
          if (grepl("flagellarpresence", current_name)) {
            meta_cell <- refined_table[current_line, "flagella"]
            refined_table[current_line, "flagella"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get growth in groups
          if (grepl("growth_in_groups", current_name)) {
            meta_cell <- refined_table[current_line, "growth in groups"]
            refined_table[current_line, "growth in groups"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get halophilic
          if (grepl("halophilic", current_name)) {
            meta_cell <- refined_table[current_line, "halophilic"]
            refined_table[current_line, "halophilic"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get pathogenic
          if ((
            !grepl("ecosystem", current_name) &
            !grepl("habitat", current_name) &
            grepl("host", current_name)
          ) |
          grepl("pathogen", current_name)) {
            meta_cell <- refined_table[current_line, "pathogenic"]
            refined_table[current_line, "pathogenic"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get metabolism
          if (grepl("metabolism", current_name)) {
            meta_cell <- refined_table[current_line, "metabolism"]
            refined_table[current_line, "metabolism"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get mobility
          if (grepl("mobility", current_name)) {
            meta_cell <- refined_table[current_line, "mobility"]
            refined_table[current_line, "mobility"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get motility
          if (grepl("motility", current_name)) {
            meta_cell <- refined_table[current_line, "motility"]
            refined_table[current_line, "motility"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get oxygen requirement
          if (grepl("oxygenreq", current_name)) {
            meta_cell <- refined_table[current_line, "oxygen requirement"]
            refined_table[current_line, "oxygen requirement"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get radioresistance
          if (grepl("radioresistance", current_name)) {
            meta_cell <- refined_table[current_line, "radioresistance"]
            refined_table[current_line, "radioresistance"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get shape
          if (grepl("shape", current_name)) {
            meta_cell <- refined_table[current_line, "shape"]
            refined_table[current_line, "shape"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get sporulation
          if (grepl("sporulation", current_name)) {
            meta_cell <- refined_table[current_line, "sporulation"]
            refined_table[current_line, "sporulation"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
          
          # get temperature range
          if (grepl("temperaturerange", current_name)) {
            meta_cell <- refined_table[current_line, "temperature range"]
            refined_table[current_line, "temperature range"] <-
              AddToCell(meta_cell,  GetTrait(names(raw_table), current_column))
          }
        }
      }
      current_column <- current_column + 1
    }
    
    # get remaining properties
    current_line <- current_line + 1
  }
  if (save_file) {
    write.csv(refined_table,
              paste0("ProTrait_v", Sys.Date(), ".csv"),
              row.names = FALSE)
  }
  message("Completed pro trait data grab")
  return(refined_table)
}
