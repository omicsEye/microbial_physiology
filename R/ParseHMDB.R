library(XML)
library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(rlist)

# configure for parallelism 
cl <- parallel::makeCluster(detectCores())
doParallel::registerDoParallel(cl)

# final formatted table
df <- as.data.frame(matrix(nrow = 1, ncol = 0))

# For each node, load list
a <- list()
# count how many nodes have been added to list
count <- 0

# For each node, flatten into one list 
flatten <- function(x) {
  m <- c()
  for (i in 1:length(x)) {
    if (!is.null(x[[i]])) {
      if (is.list(x[[i]])) {
        r <- flatten(x[[i]])
      } else {
        r <- as.list(x[[i]])
        names(r) <- names(x)[i]
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
    flatten(x[[i]])
  },
  fill = TRUE))
}

# extract necessary data from each metabolite
getMetabolite <- function(x) {
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

# Use event-driven SAX parser to process the XML without requiring the full tree structure to be loaded into memory
# Call the function defined above
start_time <- Sys.time()
xmlEventParse(
  file = 'Data/Metabolites/hmdb_metabolites.xml',
  handlers = NULL,
  trim = TRUE,
  ignoreBlanks = TRUE,
  branches = list(metabolite = getMetabolite)
)
df <- rbindlist(list(df, multi.flatten(a, flatten)), fill = TRUE)
print(Sys.time() - start_time)

write.csv(df,
          paste0("HMDB_", Sys.Date(), ".csv"),
          row.names = FALSE)

stopCluster(cl)
