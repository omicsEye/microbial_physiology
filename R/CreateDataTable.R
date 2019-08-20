#' Create microbe data table
#' Creates a single, formatted table from multiple microbial databases
#' @return a data frame
#' @export
create.microbe.data.table <- function() {
  # Source microbe information
  ##########################################
  # Protraits
  protrait <- clean.protrait()
  
  # IJSEM
  ijsem <- parse.ijsem()
  
  # BacMap
  bacmap <- bacmap.crawler()
  
  # BacDive 
  bacdive <- bacdive.crawler()
  
  ##########################################
  # Combine microbe tables
  all_table <- list(protrait, ijsem, bacmap, bacdive)
  total_table <- combine.data(all_table)
}
