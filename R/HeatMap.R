library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(colorspace)

# load and prepare abundance data from .csv to use for heat map
# returns a numerical matrix 
load.abundance.data <- function(path, column = 1) {
  if(is.na(path)) {
    path <- file.choose()
  }
  abundance_table <- read.delim(path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
  abundance_table <- abundance_table[!is.na(abundance_table[, column]), ]
  abundance_table <- abundance_table[(abundance_table[, column]) != '', ]
  abundance_table <- abundance_table[!duplicated(abundance_table[, column]), ]
  species <- abundance_table[, column]
  abundance_table <- apply(abundance_table, 2, as.numeric)
  row.names(abundance_table) <- species
  return(abundance_table)
}


# load and prepare meta data from .csv to use for annotating heat map 
# returns a data frame
load.meta.data <- function(path, tax_column = 1) {
  if(is.na(path)) {
    path <- file.choose()
  }
  data <- read.delim(path, fill=NA, stringsAsFactors = FALSE, sep=',', check.names = FALSE)
  data <- data[!is.na(data[, tax_column]), ]
  data <- data[(data[, tax_column]) != '', ]
  data <- data[!duplicated(data[, tax_column]), ]
  row.names(data) <- data[, tax_column] 
  data <- data[, -c(tax_column)]
  return(data)
}

# create a correlogram given abundance data and feature metadata
create.correlogram <- function(data, feature_meta, show = TRUE, omit = TRUE) {
  data[is.na(data)] <- 0
  
  if(omit) {
    # remove incomplete data 
    feature_meta <- na.omit(feature_meta)
    
    # Filter out data and meta data that do not have corresponding part
    feature_meta <-
      # Eliminates rows in metadata not in abundance table
      feature_meta[(rownames(feature_meta) %in% rownames(data)), , drop = FALSE]
    
    data <- # Eliminates rows in abundance table not in metadata
      data[(rownames(data) %in% rownames(feature_meta)), , drop = FALSE]
  }
  
  data <- t(data)
  data <- # Get rid of 0 correlation 
    data[, apply(data, 2, var, na.rm = TRUE) != 0]
  corMat <- cor(data, method = "spearman")
  
  return(
    pheatmap(
      corMat,
      show_colnames = FALSE,
      show_rownames = FALSE,
      treeheight_row = 0,
      treeheight_col = 0,
      fontsize_row = 6,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      annotation_row = feature_meta,
      annotation_col = feature_meta,
      cellwidth = 2,
      cellheight = 1,
      scale = "none"
    )
  )
}

# create a correlogram from ajoined data tables
multi.correlogram <- function(data_tables, sample_datas) {
  
  create.multi.correlogram <- function(a_x_b, a, b) {
    
    annotation_color_count <-
      sum(sapply(a, function(x)
        length(unique(x)))) +
      sum(sapply(b, function(x)
        length(unique(x))))
    
    newColsMeta <-
      lapply(b, function(x)
        heat.colors(length(unique(x))))
    for (i in 1:length(newColsMeta)) {
      names(newColsMeta[[i]]) <- unique(b[, i])
    }
    
    newColsMicrobe <-
      lapply(a, function(x)
        rainbow_hcl(length(unique(x))))
    for (i in 1:length(newColsMicrobe)) {
      names(newColsMicrobe[[i]]) <- unique(a[, i])
    }
    
    myColors <- append(newColsMeta, newColsMicrobe)
    
    return(
      pheatmap(
        a_x_b,
        show_colnames = FALSE,
        show_rownames = FALSE,
        treeheight_row = 0,
        treeheight_col = 0,
        fontsize_row = 6,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        annotation_row = a,
        annotation_col = b,
        # annotation_colors = myColors,
        cellwidth = 2,
        cellheight = 1,
        scale = "none"
      )
    )
  }
  
  unrefined_sample_datas <- sample_datas
  
  # check that user inputed valid data
  if (!is.list(data_tables))
    stop("Was not given list of data tables")
  if (!is.list(sample_datas))
    stop("was not given list of sample data")
  
  # create single dataframe from list of data frames for abudnance tables
  data_tables_names <- lapply(data_tables, row.names)
  data_tables <- lapply(data_tables, as.data.frame)
  data_tables <-
    as.data.frame(rbindlist(data_tables, fill = TRUE, use.names = TRUE))
  row.names(data_tables) <- unlist(data_tables_names)
  
  # create single dataframe from list of data frames for meta data
  sample_datas_names <- lapply(sample_datas, row.names)
  sample_datas <- lapply(sample_datas, as.data.frame)
  sample_datas <-
    as.data.frame(rbindlist(sample_datas, fill = TRUE, use.names = TRUE))
  row.names(sample_datas) <- unlist(sample_datas_names)
  
  # Filter out data and meta data that do not have corresponding part
  sample_datas <-
    # Eliminates rows in metadata not in abundance table
    sample_datas[(rownames(sample_datas) %in% rownames(data_tables)), , drop = FALSE]
  
  data_tables <- # Eliminates rows in abundance table not in metadata
    data_tables[(rownames(data_tables) %in% rownames(sample_datas)), , drop = FALSE]
  
  data_tables[is.na(data_tables)] <- 0
  
  data_tables <- t(data_tables)
  data_tables <- # Get rid of 0 correlation
    data_tables[, apply(data_tables, 2, var, na.rm = TRUE) != 0]
  
  corMat <- cor(data_tables, method = "spearman")
  
  # create heatmaps from all possible pairings of abundance tables 
  mat <- list()
  for(i in 1:length(sample_datas_names)) {
    for (j in 1:length(sample_datas_names)) {
      mat[[(i - 1) * length(sample_datas_names) + j]] <-
        create.multi.correlogram(corMat[row.names(corMat) %in% sample_datas_names[[i]], colnames(corMat) %in% sample_datas_names[[j]]],
                                 unrefined_sample_datas[[i]],
                                 unrefined_sample_datas[[j]])
    }
  }
  return(mat)
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
