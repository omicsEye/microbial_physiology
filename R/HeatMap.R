library(pheatmap)
library(RColorBrewer)

CreateHeatMap <- function(data, row_meta, column_meta) {
  # Filter out missing data
  row_meta <- na.omit(row_meta)
  column_meta <- na.omit(column_meta)
  data <- na.omit(data)
  
  # Filter out data and meta data that do not have corresponding part
  row_meta <- # Eliminates rows in metadata not in abundance table
    row_meta[(rownames(row_meta) %in% rownames(data)), ]
  
  data <- # Eliminates rows in abundance table not in metadata
    data[(rownames(data) %in% rownames(row_meta)),]
  
  column_meta <-
    # Eliminates rows in metadata not in abundance table
    column_meta[(rownames(column_meta) %in% colnames(data)), ]
  
  data <- data[, (colnames(data) %in% rownames(column_meta))]
  
  data <- sqrt(data)
  
  breaksList = seq(min(data) - min(data) / 500,
                   max(data) + max(data) / 500,
                   by = max(data) / 500)
  
  pheatmap(
    data.matrix(data),
    annotation_row = row_meta,
    annotation_col = column_meta,
    breaks = breaksList,
    color = colorRampPalette(rev(brewer.pal(
      n = 7, name = "RdYlBu"
    )))(length(breaksList)),
    show_colnames = FALSE,
    show_rownames = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    fontsize_row = 4,
    border = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
}
