library(XML)
library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(utils)
library(rlist)

parse.hmdb <-
  function(file = 'hmdb_metabolites.xml', link = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip') {
    
    # For each node, flatten into one list
    flatten <- function(x, pre_name = '') {
      m <- c()
      for (i in 1:length(x)) {
        if (!is.null(x[[i]])) {
          if (is.list(x[[i]])) {
            r <- flatten(x[[i]], names(x)[i])
          } else {
            r <- as.list(x[[i]])
            names(r) <- paste0(pre_name, '.', names(x)[i])
          }
          old_names <- names(m)
          m <- append(m, r)
          names(m) <- append(old_names, names(r))
        }
      }
      return(m[!duplicated(names(m))])
    }
    
    # parallel flattening
    multi.flatten <- function(x, flatten) {
      return(rbindlist(foreach(i = 1:length(x), flatten) %dopar% {
        flatten(x[[i]], '')
      },
      fill = TRUE))
    }
    
    # extract necessary data from each metabolite
    get.metabolite <- function(x) {
      a <<- list.append(a, xmlToList(x))
      remove(x)
      count <<- count + 1
      if (count >= 1000) {
        print('combining')
        df <<-
          rbindlist(list(df, multi.flatten(a, flatten)), fill = TRUE)
        a <<- list()
        count <<- 0
        
        write.csv(df,
                  paste0("HMDB_", Sys.Date(), ".csv"),
                  row.names = FALSE)
      }
    }
    
    # return number of times that pattern appears
    occurrence.symbol <- function(x, pattern) {
      return(sum(grepl(pattern, unlist(
        strsplit(x, split = NULL)
      ))))
    }
    
    # determine traits about molecular structure
    # returns a list with the respective data
    get.molecular <- function(SMILES) {
      molecular_data <- as.list(matrix(nrow = 1, ncol = 6))
      names(molecular_data) <- c(
        'number of carbons',
        'double bonds',
        'triple bonds',
        'quadruple bonds',
        'positive charge',
        'negative charge'
      )
      molecular_data['number of carbons'] <-
        occurrence.symbol(SMILES, 'C')
      molecular_data['double bonds'] <-
        occurrence.symbol(SMILES, '=')
      molecular_data['triple bonds'] <-
        occurrence.symbol(SMILES, '#')
      molecular_data['quadruple bonds'] <-
        occurrence.symbol(SMILES, '$')
      molecular_data['positive charge'] <-
        occurrence.symbol(SMILES, '+')
      molecular_data['negative charge'] <-
        occurrence.symbol(SMILES, '-')
      
      return(molecular_data)
    }
    
    # check if local file exists, download if does not
    if (!file.exists(file)) {
      download.file(link,
                    file)
    }
    
    # configure for parallelism
    cl <- parallel::makeCluster(detectCores())
    doParallel::registerDoParallel(cl)
    
    # final formatted table
    df <- as.data.frame(matrix(nrow = 1, ncol = 0))
    
    # For each node, load list
    a <- list()
    # count how many nodes have been added to list
    count <- 0
    
    # Use event-driven SAX parser to process the XML without requiring the full tree structure to be loaded into memory
    # Call the function defined above
    
    xmlEventParse(
      file = file,
      handlers = NULL,
      trim = TRUE,
      ignoreBlanks = TRUE,
      branches = list(metabolite = get.metabolite)
    )
    
    df <-
      rbindlist(list(df, multi.flatten(a, flatten)), fill = TRUE)
    
    # get molecular properties
    molecular <- rbindlist(lapply(df$'.smiles', get.molecular))
    df <- cbind(df, molecular)
    
    write.csv(df,
              paste0("HMDB_", Sys.Date(), ".csv"),
              row.names = FALSE)
    
    stopCluster(cl)
  }
