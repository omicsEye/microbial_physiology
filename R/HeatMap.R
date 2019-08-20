library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(colorspace)
library(foreach)
library(parallel)
library(doParallel)

# save a figure
save.figure <- function(figure, file_location = '') {
  ggsave(
    filename = paste0(file_location, deparse(substitute(figure)), '.png'),
    figure,
    width = 5,
    height = 6,
    units = 'in',
    dpi = 300
  )
}

# a for row, b for column
create.unique.color <- function(a, b) {
  color <-
    grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

  annotation_color_count <-
    sum(sapply(a, function(x)
      length(unique(x)))) +
    sum(sapply(b, function(x)
      length(unique(x))))

  newColsMeta <-
    lapply(b, function(x)
      color[seq(18, length(color), floor(length(color) / length(na.omit(unique(x)))))])
  for (i in 1:length(newColsMeta)) {
    names(newColsMeta[[i]]) <- unique(na.omit(b[, i]))
  }

  newColsMicrobe <-
    lapply(a, function(x)
      color[seq(18, length(color), floor(length(color) / length(na.omit(unique(x)))))])
  for (i in 1:length(newColsMicrobe)) {
    names(newColsMicrobe[[i]]) <- unique(na.omit(a[, i]))
  }
  return(append(newColsMeta, newColsMicrobe))
}

# create a correlogram given abundance data and feature metadata
create.correlogram <-
  function(data,
           feature_meta,
           show = TRUE,
           omit = TRUE,
           cluster_distance_method = 'euclidean') {
    data[is.na(data)] <- 0

    if (omit) {
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
        clustering_distance_rows = cluster_distance_method,
        clustering_distance_cols = cluster_distance_method,
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        annotation_row = feature_meta,
        annotation_col = feature_meta,
        annotation_colors = create.unique.color(feature_meta, feature_meta),
        cellwidth = 2,
        cellheight = 1,
        scale = "none",
        silent = !show,
        legend = FALSE
      )
    )
  }

# create a correlogram from ajoined data tables
multi.correlogram <- function(data_tables, sample_datas) {
  create.multi.correlogram <- function(a_x_b, a, b) {
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

  data_tables <-
    # Eliminates rows in abundance table not in metadata
    data_tables[(rownames(data_tables) %in% rownames(sample_datas)), , drop = FALSE]

  data_tables[is.na(data_tables)] <- 0

  data_tables <- t(data_tables)
  data_tables <- # Get rid of 0 correlation
    data_tables[, apply(data_tables, 2, var, na.rm = TRUE) != 0]

  corMat <- cor(data_tables, method = "spearman")

  # create heatmaps from all possible pairings of abundance tables
  mat <- list()
  for (i in 1:length(sample_datas_names)) {
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
  function(data,
           sample_meta,
           feature_meta,
           percentile = 0.75,
           filter = '',
           show = FALSE,
           omit_na = TRUE,
           cluster_distance_method = 'euclidean') {
    # Filter out missing data
    data[is.na(data)] <- 0

    if (omit_na) {
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
    }

    if (filter == 'abundance') {
      data <-
        data[rowMeans(data) >= quantile(rowMeans(data), percentile),]
    } else if (filter == 'variance') {
      data <-
        data[apply(data, 1, var) >= quantile(apply(data, 1, var), percentile),]
    }

    # data[data > 0] <- sqrt(data[data > 0])
    # data[data < 0] <- sqrt(abs(data[data < 0])) * -1
    data <- data / max(data)

    breaksList = seq(min(data) - min(data) / 500,
                     max(data) + max(data) / 500,
                     by = max(data) / 500)

    return(
      pheatmap(
        data,
        annotation_row = feature_meta,
        annotation_col = sample_meta,
        annotation_colors = create.unique.color(feature_meta, sample_meta),
        breaks = breaksList,
        color = colorRampPalette(rev(brewer.pal(
          n = 7, name = "RdYlBu"
        )))(length(breaksList)),
        show_colnames = FALSE,
        show_rownames = FALSE,
        treeheight_row = 0,
        treeheight_col = 0,
        fontsize = 4,
        border = TRUE,
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        clustering_distance_rows = cluster_distance_method,
        clustering_distance_cols = cluster_distance_method,
        silent = !show
      )
    )
  }

# create a heatmap for each morphological category
all.heatmap <-
  function(data,
           sample_meta,
           feature_meta,
           percentile = 0.75,
           show = FALSE,
           omit_na = FALSE,
           file_location = '',
           cluster_distance_method = "euclidean") {
    foreach(i = 1:length(feature_meta)) %do% {
      ggsave(
        filename = paste0(file_location, names(feature_meta)[i], '.png'),
        create.heatmap(
          data = data,
          sample_meta = sample_meta,
          feature_meta = feature_meta[, i, drop = FALSE],
          percentile = percentile,
          show = show,
          omit_na = omit_na,
          cluster_distance_method = cluster_distance_method
        ),
        width = 5,
        height = 6,
        units = 'in',
        dpi = 300
      )
    }
  }

# create a heat map of a single feature
one.v.all <-
  function(data,
           sample_meta,
           feature_meta,
           which = 2,
           percentile = 0.75,
           show = FALSE,
           column,
           trait,
           cluster_distance_method = "euclidean") {
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
    return(
      create.heatmap(
        data,
        sample_meta,
        feature_meta,
        percentile = percentile,
        show = show,
        omit_na = FALSE,
        cluster_distance_method = cluster_distance_method
      )
    )
  }

# create all possible heatmaps for a certain feature
all.one.v.all <-
  function(data,
           sample_meta,
           feature_meta,
           which = 2,
           percentile = 0.75,
           show = FALSE,
           column,
           directory = '',
           cluster_distance_method = "euclidean") {
    if (which == 1) {
      all_observers <-
        unique(unlist(strsplit(as.character(sample_meta[, column]), ",")))
    } else {
      all_observers <-
        unique(unlist(strsplit(as.character(feature_meta[, column]), ",")))
    }
    all_observers <- all_observers[!is.na(all_observers)]

    create <- function(x) {
      plot <-
        one.v.all(
          data,
          sample_meta,
          feature_meta,
          percentile = percentile,
          show = show,
          which = which,
          column = column,
          trait = x,
          cluster_distance_method = cluster_distance_method
        )
      ggsave(
        file = paste0(directory, x, '_v_all_', Sys.Date(), '.png'),
        plot = plot,
        width = 5,
        height = 4,
        units = 'in',
        dpi = 300
      )
    }

    lapply(all_observers, create)
  }
