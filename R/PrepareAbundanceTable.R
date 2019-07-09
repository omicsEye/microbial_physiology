
prepare.abundance.table <- function(file, tax = 's', tax_col = 1) {
  # available taxanomic information
  available_tax <- c('k', 'p', 'c', 'o', 'f', 'g', 's', 't')
  if (!(tax %in% available_tax))
    stop(paste0(
      'invalid taxonomic abbreviation given, available are: ',
      toString(available_tax, collapse = ', ')
    ))
  
  # determine taxonomic subcategory below user specified
  next_tax <- match(tax, available_tax) + 1
  if (next_tax <= length(available_tax))
    next_tax <- paste0(available_tax[next_tax], '__')
  tax <- paste0(tax, '__')
  
  # load abundance table
  if (!file.exists(file))
    stop(paste0('File: ', file, " not found"))
  abundance_table <-
    read.table(
      file,
      header = TRUE,
      comment.char = "",
      check.names = FALSE
    )
  
  if (tax_col < 0 | tax_col > ncol(abundance_table))
    stop('Taxonomy column given not within the dimension of the abundance table')
  
  # for each in specie column, eliminate if not designated type
  # Determine if want kingdom, phylum, class, order, family, genus, or specie (default)
  abundance_table <-
    abundance_table[grepl(tax, abundance_table[, tax_col]), ]
  
  # delete all that have more taxanomic information
  abundance_table <-
    abundance_table[!grepl(next_tax, abundance_table[, tax_col]), ]
  
  # process names
  process.names <- function(name) {
    name <- substring(name, regexpr(tax, name) + 3)
    return(gsub("_", " ", name))
  }
  
  process_names <- lapply(abundance_table[, tax_col], process.names)
  row.names(abundance_table) <- process_names
  
  abundance_table <- abundance_table[, -1]
  
  write.csv(abundance_table,
            paste0(file, '_AbundanceTable_', Sys.Date(), '.csv'))
  
  return(abundance_table)
}
