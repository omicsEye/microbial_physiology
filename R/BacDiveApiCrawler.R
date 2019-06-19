library(RCurl)
library(rjson)

# get JSON response from a given url
GetJSONResponse <- function(usrname, pass, url) {
  # get server response
  response <-
    getURL(url,
           userpwd = paste0(usrname, ":", pass),
           httpauth = 1L)
  
  # convert to JSON format
  jsondata <- fromJSON(response)
  
  return(jsondata)
}

# given a cell entry, update it based on a new entry
UpdateCell <- function(current_entry, new_entry) {
  if (length(new_entry) > 0) {
    # update if no previous entry
    if (length(current_entry) <= 0 || is.na(current_entry)) {
      return(new_entry)
      
      # update if new entry is novel
    } else if (!grepl(new_entry, current_entry)) {
      return(paste(current_entry, new_entry, sep = ', '))
      
      # do not update if new entry not novel
    } else {
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

# update a cell given data that is expected to have multiple entries
GetMultiData <- function(current_entry, data, data_entry, checker) {
  # parse through all entries
  for (j in 1:length(data)) {
    current_sample <- data[[j]]
    ability <-
      current_sample[checker]  # determines whether pocesses trait
    
    # add trait if known or unknown to pocess
    if (length(ability) > 0 &&
        (ability == 'positive' ||
         ability == "+" || ability == 'TRUE') ||
        length(ability) <= 0) {
      new_entry <- current_sample[data_entry]
      current_entry <- UpdateCell(current_entry, new_entry)
    }
  }
  return(current_entry)
}

# Queries BacDove API
BacDiveCrawler <- function(usrname, pass, save_file = TRUE) {
  # All traits that are to be extracted from BacDive data
  taxonomic_trait <-
    c(
      "Species",
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
      "Flagella",
      "Incubation period",
      "Ability of spore formation",
      "Type of spore",
      "Multicellular complex forming ability",
      "Murein short key",
      "Murein types",
      "Oxygen tolerance",
      "Nutrition type",
      "met_antibiotica",
      "halophily",
      "met_util",
      "met_production",
      "enzymes",
      "Temperature range",
      "pH",
      "Pathogenic in humans",
      "Pathogenic in animals",
      "Pathogenic in plants"
    )
  
  # table that will hold all extracted data
  table <-
    data.frame(matrix(NA, nrow = 0, ncol = length(taxonomic_trait)))
  names(table) <- taxonomic_trait
  
  # to find links to microbial data, this is the first such page
  next_page <-
    'https://bacdive.dsmz.de/api/bacdive/bacdive_id/?page=1'
  
  # continue running until there is no page left
  while (length(next_page) > 0 || !is.na(next_page)) {
    # get all urls from page
    allbacdiveIDs <- next_page
    url_allbacdiveIDs <-
      URLencode(paste0(allbacdiveIDs, '&format=json'))
    
    id_result <- GetJSONResponse(usrname, pass, url_allbacdiveIDs)
    next_page <-
      id_result$'next'  # update what the next search page will be
    results <- id_result$results  # all microbial links on the page
    
    # querry all microbial links found
    for (i in 1:length(results)) {
      # create a row for the new microbe
      current_row <-
        data.frame(matrix(NA, nrow = 1, ncol = length(taxonomic_trait)))
      names(current_row) <- taxonomic_trait
      
      # get the data about the microbe
      current_search <- results[[i]][[1]]
      url_bacdiveID <-
        URLencode(paste0(current_search, '/?format=json'))
      species_result <-
        GetJSONResponse(usrname, pass, url_bacdiveID)
      
      # get taxonomy data
      taxonomy_data <- species_result$taxonomy_name$strains[[1]]
      current_row$Domain <- UpdateCell(NA, taxonomy_data$domain)
      current_row$Phylum <- UpdateCell(NA, taxonomy_data$phylum)
      current_row$Class <- UpdateCell(NA, taxonomy_data$class)
      current_row$Order <- UpdateCell(NA, taxonomy_data$ordo)
      current_row$Family <- UpdateCell(NA, taxonomy_data$family)
      current_row$Genus <- UpdateCell(NA, taxonomy_data$genus)
      current_row$Species <- UpdateCell(NA, taxonomy_data$species)
      
      # get morhology data
      morphology_data <- species_result$morphology_physiology
      
      # cell morphology
      cell_data <- morphology_data$cell_morphology[[1]]
      current_row$`Gram stain` <-
        UpdateCell(NA, cell_data$gram_stain)
      current_row$`Cell length` <-
        UpdateCell(NA, cell_data$cell_len)
      current_row$`Cell width` <-
        UpdateCell(NA, cell_data$cell_width)
      current_row$`Cell shape` <-
        UpdateCell(NA, cell_data$cell_shape)
      current_row$Motility <- UpdateCell(NA, cell_data$motility)
      current_row$Flagella <-
        UpdateCell(NA, cell_data$flagellum_arrangement)
      
      # colony data
      colony_data <- morphology_data$colony_morphology[[1]]
      current_row$`Type of hemolysis` <-
        UpdateCell(NA, colony_data$hemolysis_type)
      current_row$`Hemolysis Ability` <-
        UpdateCell(NA, colony_data$hemolysis_ability)
      current_row$`Colony size` <-
        UpdateCell(NA, colony_data$colony_len)
      current_row$`Colony color` <-
        UpdateCell(NA, colony_data$colony_color)
      current_row$`Colony shape` <-
        UpdateCell(NA, colony_data$colony_shape)
      current_row$`Incubation period` <-
        UpdateCell(NA, colony_data$incubation_period)
      
      # mutlicellular data
      multicellular_data <-
        morphology_data$multicellular_morphology[[1]]
      current_row$`Multicellular complex forming ability` <-
        UpdateCell(NA, multicellular_data$ability)
      
      # spore data
      spore_data <- morphology_data$spore_formation[[1]]
      current_row$`Type of spore` <- UpdateCell(NA, spore_data$type)
      current_row$`Ability of spore formation` <-
        UpdateCell(NA, spore_data$ability)
      
      # murein data
      murein_data <- morphology_data$murein[[1]]
      current_row$`Murein short key` <-
        UpdateCell(NA, murein_data$murein_short_index)
      current_row$`Murein types` <-
        UpdateCell(NA, murein_data$murein_types)
      
      # nutrition data
      nutrition_data <- morphology_data$nutrition_type[[1]]
      current_row$`Nutrition type` <-
        UpdateCell(NA, nutrition_data$nutrition_type)
      
      # oxygen tolerance data
      oxygen_data <- morphology_data$oxygen_tolerance[[1]]
      current_row$`Oxygen tolerance` <-
        UpdateCell(NA, oxygen_data$oxygen_tol)
      
      # compound production data
      compound_data <- morphology_data$compound_production
      current_row$met_production <-
        GetMultiData(current_row$met_production,
                     compound_data,
                     'compound_name',
                     NULL)
      compound_data <- morphology_data$met_production
      current_row$met_production <-
        GetMultiData(current_row$met_production,
                     compound_data,
                     'metabolite_prod',
                     'production')
      
      # metabolite utilization data
      compound_data <- morphology_data$met_util
      current_row$met_util <-
        GetMultiData(current_row$met_util,
                     compound_data,
                     'metabolite_util',
                     'ability')
      
      # antiobiotic resistance data
      anti_data <- morphology_data$met_antibiotica
      current_row$met_antibiotica <-
        GetMultiData(current_row$met_antibiotica,
                     anti_data,
                     'metabolite_antib',
                     'ab_resistant')
      
      
      # enzyme data
      enzyme_data <- morphology_data$enzymes
      current_row$enzymes <-
        GetMultiData(current_row$enzymes, enzyme_data, 'enzyme', 'activity')
      
      # halophily data
      halophily_data <- morphology_data$halophily
      current_row$halophily <-
        GetMultiData(current_row$halophily,
                     halophily_data,
                     "salt_concentration",
                     'ability')
      
      # get temperature data
      temperature_data <-
        species_result$culture_growth_condition$culture_temp
      current_row$`Temperature range` <-
        GetMultiData(
          current_row$`Temperature range`,
          temperature_data,
          'temperature_range',
          'ability'
        )
      
      # get pH data
      ph_data <- species_result$culture_growth_condition$culture_pH
      current_row$pH <-
        GetMultiData(current_row$pH, ph_data, 'pH', 'ability')
      
      # get pathogenic data
      pathogenic_data <-
        species_result$application_interaction$risk_assessment[[1]]
      current_row$`Pathogenic in humans` <-
        UpdateCell(NA, pathogenic_data$pathogenicity_human)
      current_row$`Pathogenic in animals` <-
        UpdateCell(NA, pathogenic_data$pathogenicity_animal)
      current_row$`Pathogenic in plants` <-
        UpdateCell(NA, pathogenic_data$pathogenicity_plant)
      
      # update tabe with new row
      table <- rbind(table, current_row)
    }
    message("Gathered 100 additional results")  # prompt user about progress
  }
  table <-
    apply(table, 2, as.character)  # remove unwanted formating
  
  # save table to working directory
  if (save_file) {
    write.csv(table, paste0("BacDive_v", Sys.Date(), ".csv"), row.names = FALSE)
  }
  message("Completed bacdive query")  # prompt user about completion
  return(table)
}
