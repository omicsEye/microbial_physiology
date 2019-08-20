#' Enrichment Analysis
#' Determines ratio of trait possessing to non-trait possessing across all samples
#' 
#' @param data_table abundance table of samples
#' @param meta_data table of types and their traits
#' @param category which trait type to perform enrichment on
#' @param trait specific trait to determine if possesses
#' @param na_omit whether to eliminate all species that do not know trait of
#' 
#' @return ratios 
#' @export
enrichmentAnalysis <- function(data_table, meta_data, category, trait, na_omit = TRUE) {
  
  meta_data <- meta_data[, category, drop = FALSE]
  # remove unlabled microbes
  if (na_omit)
    meta_data <- na.omit(meta_data)
  data_table <- data_table[row.names(data_table) %in% row.names(meta_data), ]
  # removed unneeded microbes
  meta_data <- meta_data[row.names(meta_data) %in% row.names(data_table), , drop = FALSE]
  # order metadata same as abundance table
  meta_data <- meta_data[order(match(rownames(meta_data), rownames(data_table))), , drop = FALSE]
  
  ratios <- c(NULL)
  for (i in 1:ncol(data_table)) {
    # join metadata with sample
    current_sample <- cbind(data_table[, i, drop = FALSE], meta_data, stringsAsFactors = FALSE)
    
    # remove all microbes that do not exist in sample
    current_sample <- current_sample[which(current_sample[, 1] != 0), , drop = FALSE]
    
    if (nrow(current_sample) > 0 ) {
      has_sum <- 0
      has_not_sum <- 0
      for (j in 1:nrow(current_sample)) {
        if (!is.na(current_sample[j, category]) && current_sample[j, category] == trait) {
          has_sum <- has_sum + current_sample[j, 1]
        } else {
          has_not_sum <- has_not_sum + current_sample[j, 1]
        }
      }
      ratios <- c(ratios, (has_sum / has_not_sum))
    }
  }
  return(ratios)
}


# densityOf <- density(ratios)
# plot(densityOf, xlim = c(min(densityOf$x), 10))
# abline(v = 1, lty = 3)
# pValue <- t.test(ratios, mu = 1)
