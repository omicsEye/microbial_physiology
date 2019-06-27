library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# create a correlogram given abundance data and feature metadata
create.correlogram <- function(data, feature_meta, show = TRUE) {
  # remove incomplete data 
  feature_meta <- na.omit(feature_meta)
  data <- na.omit(data)
  
  # Filter out data and meta data that do not have corresponding part
  feature_meta <-
    # Eliminates rows in metadata not in abundance table
    feature_meta[(rownames(feature_meta) %in% rownames(data)), , drop = FALSE]
  
  data <- # Eliminates rows in abundance table not in metadata
    data[(rownames(data) %in% rownames(feature_meta)), , drop = FALSE]
  
  data <- t(data)
  data <- # Get rid of 0 correlation 
    data[, apply(data, 2, var, na.rm = TRUE) != 0]
  corMat <- cor(data, method = "spearman")
  
  breaksList = seq(min(data) - min(data) / 500,
                   max(data) + max(data) / 500,
                   by = max(data) / 500)
  
  return(
    pheatmap(
      data.matrix(corMat),
      annotation_row = feature_meta,
      annotation_col = feature_meta,
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
      clustering_distance_cols = "euclidean",
      silent = !show
    )
  )
}

# create a heat map given abundance data, sample metadata, and feature metadata
create.heatmap <-
  function(data, sample_meta, feature_meta, percentile = 0.75, show = FALSE, omit_na = TRUE) {
    # Filter out missing data
    if (omit_na) {
      feature_meta <- na.omit(feature_meta)
      sample_meta <- na.omit(sample_meta)
      data <- na.omit(data)
    }
    
    data <-
      data[rowMeans(data) >= quantile(rowMeans(data), percentile), ] # Filtering for percentile
    
    # Filter out data and meta data that do not have corresponding part
    feature_meta <-
      # Eliminates rows in metadata not in abundance table
      feature_meta[(rownames(feature_meta) %in% rownames(data)), , drop = FALSE]
    
    data <- # Eliminates rows in abundance table not in metadata
      data[(rownames(data) %in% rownames(feature_meta)), , drop = FALSE]
    
    sample_meta <-
      # Eliminates rows in metadata not in abundance table
      sample_meta[(rownames(sample_meta) %in% colnames(data)),  , drop = FALSE]
    
    data <-
      data[, (colnames(data) %in% rownames(sample_meta)), drop = FALSE]
    
    data <- sqrt(data)
    
    breaksList = seq(min(data) - min(data) / 500,
                     max(data) + max(data) / 500,
                     by = max(data) / 500)
    
    return(
      pheatmap(
        data.matrix(data),
        annotation_row = feature_meta,
        annotation_col = sample_meta,
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
        clustering_distance_cols = "euclidean",
        silent = !show
      )
    )
  }

# create a heat map of a single feature
one.v.all <-
  function(data, sample_meta, feature_meta, which = 2, percentile = 0.75, show = FALSE, column, trait) {
    meta <- ""
    
    if (which == 1) {
      meta <-
        as.data.frame(sample_meta[, column], row.names = rownames(sample_meta))
    } else {
      meta <-
        as.data.frame(feature_meta[, column], row.names = rownames(feature_meta))
    }
    names(meta) <- column
    
    get.trait <- function(x) {
      if (!grepl(trait, x, fixed = TRUE)) {
        return('other')
      } else {
        return(trait)
      }
    }
    
    raw_meta <- lapply(meta[, column], get.trait)
    
    meta[, column] <- t(as.data.frame(raw_meta))
    
    
    if (which == 1) {
      sample_meta <- meta
    } else {
      feature_meta <- meta
    }
    return(create.heatmap(data, sample_meta, feature_meta, percentile = percentile, show = show, omit_na = FALSE))
  }

# create all possible heatmaps for a certain feature
all.one.v.all <- 
  function(data, sample_meta, feature_meta, which = 2, percentile = 0.75, show = FALSE, column, directory='') {
    
    all_observers <- 
      unique(unlist(strsplit(as.character(feature_meta[, column]), ",")))
    all_observers <- all_observers[!is.na(all_observers)]
    
    create <- function(x) {
      plot <- one.v.all(data, sample_meta, feature_meta, percentile = percentile, show = show, which = which, column, x)
      ggsave(
        file = paste0(directory, x, '_v_all_', Sys.Date(), '.png'), 
        plot = plot,
        width = 5,
        height = 4,
        units = 'in',
        dpi = 300)
    }
    
    lapply(all_observers, create)
  }
