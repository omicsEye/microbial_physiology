library(rjson)
library(R.utils)
library(crul)
library(RCurl)
library(jsonlite)

source('R/Utility.R')

# get JSON response from a given url
get.json.response <- function(usrname, pass, url) {
  # get server response
  response <-
    getURL(url,
           userpwd = paste0(usrname, ":", pass),
           httpauth = 1L)
  
  # convert to JSON format
  jsondata <- fromJSON(response)
  
  return(jsondata)
}

# update a cell given data that is expected to have multiple entries
get.multi.data <-
  function(current_entry, data, data_entry, checker) {
    if (is.data.frame(data)) {
      # parse through all entries
      for (j in 1:nrow(data)) {
        current_sample <- data[j, , drop = FALSE]
        ability <-
          current_sample[checker]  # determines whether pocesses trait
        
        # add trait if known or unknown to pocess
        if (length(ability) > 0 && !is.na(ability) &&
            (ability == 'positive' ||
             ability == "+" || ability == 'TRUE') ||
            length(ability) <= 0) {
          new_entry <- current_sample[data_entry]
          current_entry <- update.cell(current_entry, new_entry)
        }
      }
    }
    return(current_entry)
  }

# Queries BacDove API
bacdive.crawler <-
  function(usrname = 'mbarbini@broadinstitute.org',
           pass = 'trellointeractions',
           num_requests = 10,
           save_file = TRUE) {
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
        "Antibiotic Resistance",
        "Antibiotic Sensitivity",
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
    
    # table that will hold all extracted data
    table <-
      data.frame(matrix(NA, nrow = 0, ncol = length(taxonomic_trait)))
    names(table) <- taxonomic_trait
    
    # to determine how many microbes should find
    count_page <-
      'https://bacdive.dsmz.de/api/bacdive/bacdive_id/?page=1/&format=json'
    response <- get.json.response(usrname, pass, count_page)
    count <- response$count
    
    # to find links to microbial data, this is the base url
    url_page <-
      'https://bacdive.dsmz.de/api/bacdive/bacdive_id/'
    
    page <- 1  # which microbe page currrently on
    
    message("Began BacDive")
    
    num_results <- 0  # how many species retrieved
    num_no_find <- 0 # how many empty entries in a row
    
    while (num_results <= count) {
      url_list <- c()  # which urls to request
      
      for (i in 1:num_requests) {
        url_list <-
          append(url_list, URLencode(paste0(url_page, page, '/?format=json')))
        page <- page + 1
        
        # check progress
        if (page %% 100 == 0) {
          num_results <- nrow(table)
          message(paste("Gathered", num_results, "results"))  # prompt user about progress
          
          # save copy in case of crash
          if (save_file &&
              num_results %% 5000 < 100 && num_results > 0) {
            table <-
              apply(table, 2, as.character)  # remove unwanted formating
            
            write.csv(table,
                      paste0("BacDive_v", Sys.Date(), ".csv"),
                      row.names = FALSE)
          }
        }
      }
      
      # download URLs
      responses <-
        try(withTimeout(
          Async$new(urls = url_list,
                    auth = auth(user = usrname, pwd = pass)),
          timeout = 5,
          onTimeout = "error"
        ))
      responses <- responses$get()
      
      if (class(responses) != "try-error") {
        for (i in 1:length(responses)) {
          # get the data about the microbe
          
          species_result <- responses[[i]]
          if (!is.atomic(species_result)) {
            species_result <- responses[[i]]$parse("UTF-8")
            if (is.character(species_result) &&
                validate(species_result)) {
              species_result <- fromJSON(species_result)
            }
          }
          
          if (is.list(species_result) &&
              is.null(species_result$detail))  {
            # create a row for the new microbe
            current_row <-
              data.frame(matrix(
                NA,
                nrow = 1,
                ncol = length(taxonomic_trait)
              ))
            names(current_row) <- taxonomic_trait
            
            # get taxonomy data
            taxonomy_data <-
              species_result$taxonomy_name$strains
            current_row$Domain <-
              update.cell(NA, taxonomy_data$domain)
            current_row$Phylum <-
              update.cell(NA, taxonomy_data$phylum)
            current_row$Class <-
              update.cell(NA, taxonomy_data$class)
            current_row$Order <- update.cell(NA, taxonomy_data$ordo)
            current_row$Family <-
              update.cell(NA, taxonomy_data$family)
            current_row$Genus <-
              update.cell(NA, taxonomy_data$genus)
            current_row$Species <-
              update.cell(NA, taxonomy_data$species)
            
            # get morhology data
            morphology_data <- species_result$morphology_physiology
            
            # cell morphology
            cell_data <- morphology_data$cell_morphology
            current_row$`Gram stain` <-
              update.cell(NA, cell_data$gram_stain)
            current_row$`Cell length` <-
              update.cell(NA, cell_data$cell_len)
            current_row$`Cell width` <-
              update.cell(NA, cell_data$cell_width)
            current_row$`Cell shape` <-
              update.cell(NA, cell_data$cell_shape)
            current_row$Motility <-
              update.cell(NA, cell_data$motility)
            current_row$Flagella <-
              update.cell(NA, cell_data$flagellum_arrangement)
            
            # colony data
            colony_data <- morphology_data$colony_morphology
            current_row$`Type of hemolysis` <-
              update.cell(NA, colony_data$hemolysis_type)
            current_row$`Hemolysis Ability` <-
              update.cell(NA, colony_data$hemolysis_ability)
            current_row$`Colony size` <-
              update.cell(NA, colony_data$colony_len)
            current_row$`Colony color` <-
              update.cell(NA, colony_data$colony_color)
            current_row$`Colony shape` <-
              update.cell(NA, colony_data$colony_shape)
            current_row$`Incubation period` <-
              update.cell(NA, colony_data$incubation_period)
            
            # mutlicellular data
            multicellular_data <-
              morphology_data$multicellular_morphology
            current_row$`Multicellular complex forming ability` <-
              update.cell(NA, multicellular_data$ability)
            
            # spore data
            spore_data <- morphology_data$spore_formation
            current_row$`Type of spore` <-
              update.cell(NA, spore_data$type)
            current_row$`Ability of spore formation` <-
              update.cell(NA, spore_data$ability)
            
            # murein data
            murein_data <- morphology_data$murein
            current_row$`Murein short key` <-
              update.cell(NA, murein_data$murein_short_index)
            current_row$`Murein types` <-
              update.cell(NA, murein_data$murein_types)
            
            # nutrition data
            nutrition_data <- morphology_data$nutrition_type
            current_row$`Nutrition type` <-
              update.cell(NA, nutrition_data$nutrition_type)
            
            # oxygen tolerance data
            oxygen_data <- morphology_data$oxygen_tolerance
            current_row$`Oxygen tolerance` <-
              update.cell(NA, oxygen_data$oxygen_tol)
            
            # compound production data
            compound_data <- morphology_data$compound_production
            current_row$`Metabolite Production` <-
              get.multi.data(
                current_row$`Metabolite Production`,
                compound_data,
                'compound_name',
                NULL
              )
            compound_data <- morphology_data$met_production
            current_row$`Metabolite Production` <-
              get.multi.data(
                current_row$`Metabolite Production`,
                compound_data,
                'metabolite_prod',
                'production'
              )
            
            # metabolite utilization data
            compound_data <- morphology_data$met_util
            current_row$`Metabolite Utilization` <-
              get.multi.data(
                current_row$`Metabolite Utilization`,
                compound_data,
                'metabolite_util',
                'ability'
              )
            
            # antiobiotic resistance data
            anti_data <- morphology_data$met_antibiotica
            current_row$`Antibiotic Resistance` <-
              get.multi.data(
                current_row$`Antibiotic Resistance`,
                anti_data,
                'metabolite_antib',
                'ab_resistant'
              )
            
            # antibiotic sensitivity data
            current_row$`Antibiotic Sensitivity` <-
              get.multi.data(
                current_row$`Antibiotic Sensitivity`,
                anti_data,
                'metabolite_antib',
                'ab_sensitive'
              )
            
            
            # enzyme data
            enzyme_data <- morphology_data$enzymes
            current_row$enzymes <-
              get.multi.data(current_row$enzymes,
                             enzyme_data,
                             'enzyme',
                             'activity')
            
            # halophily data
            halophily_data <- morphology_data$halophily
            current_row$halophily <-
              get.multi.data(current_row$halophily,
                             halophily_data,
                             "salt_concentration",
                             'ability')
            
            # get temperature data
            temperature_data <-
              species_result$culture_growth_condition$culture_temp
            current_row$`Temperature range` <-
              get.multi.data(
                current_row$`Temperature range`,
                temperature_data,
                'temperature_range',
                'ability'
              )
            
            # get pH data
            ph_data <-
              species_result$culture_growth_condition$culture_pH
            current_row$pH <-
              get.multi.data(current_row$pH, ph_data, 'pH', 'ability')
            
            # get pathogenic data
            pathogenic_data <-
              species_result$application_interaction$risk_assessment
            current_row$`Pathogenic in humans` <-
              update.cell(NA, pathogenic_data$pathogenicity_human)
            current_row$`Pathogenic in animals` <-
              update.cell(NA, pathogenic_data$pathogenicity_animal)
            current_row$`Pathogenic in plants` <-
              update.cell(NA, pathogenic_data$pathogenicity_plant)
            
            # update table with new row
            table <- rbind(table, current_row)
            
            num_no_find <- 0
            
          } else {
            print("didnt find anything on current page")
            
            num_no_find <- num_no_find + 1
            
            if (num_no_find > 1000 && num_results > 100) {
              print("Attempting to find correct entry")
              
              all_page <- as.integer(num_results / 100)
              allbacdiveIDs <-
                'https://bacdive.dsmz.de/api/bacdive/bacdive_id/?page='
              url_allbacdiveIDs <-
                URLencode(paste0(allbacdiveIDs, all_page, '&format=json'))
              
              id_result <-
                get.json.response(url = url_allbacdiveIDs,
                                  usrname = 'mbarbini@broadinstitute.org',
                                  pass = 'trellointeractions')
              
              result <-
                id_result$results$url[[(num_results %% 100) + 1]]
              
              page <-
                as.numeric(gsub("[^\\d]+", "", result, perl = TRUE))
              num_no_find <- 0
            }
            
          }
        }
      } else {
        print("skipped URL due to timeout")
      }
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