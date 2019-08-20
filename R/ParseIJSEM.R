#' Parse IJSEM datatable
#' @details downloads a local copy of the IJSEM datatable and cleans
#'
#' @return a data frame representation of the IJSEM datatable
#' @export
parse.ijsem <- function() {
  # download data table
  if (!file.exists('IJSEM_pheno_db_v1.0.txt')) {
    download.file(
      'https://ndownloader.figshare.com/files/6994457',
      'IJSEM_pheno_db_v1.0.txt'
    )
  }
  raw_table <-
    read.delim(
      'IJSEM_pheno_db_v1.0.txt',
      fill = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  for (i in 1:nrow(raw_table)) {
    # merge species and genus name
    raw_table$`species name`[i] <-
      paste0(raw_table$`Genus name`[i], " ", raw_table$`species name`[i])
    # merge other and eco system
    if (!is.na(raw_table$Habitat[i]) && raw_table$Habitat[i] == 'other') {
      raw_table$Habitat[i] <- raw_table$`If 'other' was chosen above, please enter a habitat below`[i]
    }
  }

  write.csv(raw_table,
            paste0("IJSEM_v", Sys.Date(), ".csv"),
            row.names = FALSE)

  return(raw_table)
}
