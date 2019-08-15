

prepare.blah <-
  function(ihmp,
           microbe_abundance_table,
           metabolite_data,
           microbe_data) {
    metabolite_names <- row.names(ihmp)
    microbe_names <- row.names(microbe_abundance_table)
    
    data_tables = list(ihmp, microbe_abundance_table)
    sample_datas = list(metabolite_data, microbe_data)
    
    data_tables_names <- unlist(lapply(data_tables, row.names))
    data_tables <- lapply(data_tables, as.data.frame)
    data_tables <-
      as.data.frame(rbindlist(data_tables, fill = TRUE, use.names = TRUE))
    row.names(data_tables) <- data_tables_names
    
    sample_datas_names <- unlist(lapply(sample_datas, row.names))
    sample_datas <- lapply(sample_datas, as.data.frame)
    sample_datas <-
      as.data.frame(rbindlist(sample_datas, fill = TRUE, use.names = TRUE))
    row.names(sample_datas) <- sample_datas_names
    
    data <- data_tables
    feature_meta <- sample_datas
    show <- TRUE
    
    # Filter out data and meta data that do not have corresponding part
    feature_meta <-
      # Eliminates rows in metadata not in abundance table
      feature_meta[(rownames(feature_meta) %in% rownames(data)), , drop = FALSE]
    
    data <- # Eliminates rows in abundance table not in metadata
      data[(rownames(data) %in% rownames(feature_meta)), , drop = FALSE]
    
    data[is.na(data)] <- 0
    
    data <- t(data)
    data <- # Get rid of 0 correlation
      data[, apply(data, 2, var, na.rm = TRUE) != 0]
    corMat <- cor(data, method = "spearman")
    
    metabolite_x_metabolite <-
      corMat[row.names(corMat) %in% metabolite_names, colnames(corMat) %in% metabolite_names]
    
    microbe_x_microbe <-
      corMat[row.names(corMat) %in% microbe_names, colnames(corMat) %in% microbe_names]
    
    microbe_x_metabolite <-
      corMat[row.names(corMat) %in% microbe_names, colnames(corMat) %in% metabolite_names]
    
    return(microbe_x_metabolite)
  }

create.blah <-
  function(microbe_x_metabolite,
           microbe_meta,
           metabolite_meta) {
    
    # microbe_x_metabolite <-
    #   microbe_x_metabolite[rowMeans(microbe_x_metabolite) >= quantile(rowMeans(microbe_x_metabolite), 0.75), ] # Filtering for percentile
    # 
    breaksList = seq(min(microbe_x_metabolite) - min(microbe_x_metabolite) / 500,
                     max(microbe_x_metabolite) + max(microbe_x_metabolite) / 500,
                     by = max(abs(microbe_x_metabolite)) / 500)
    
    annotation_color_count <-
      sum(sapply(microbe_meta, function(x) length(unique(x)))) + 
      sum(sapply(metabolite_meta, function(x) length(unique(x))))
    
    newColsMeta <- lapply(metabolite_meta, function(x) heat.colors(length(unique(x))))
    for(i in 1:length(newColsMeta)) {
      names(newColsMeta[[i]]) <- unique(metabolite_meta[, i])
    }
    
    newColsMicrobe <- lapply(microbe_meta, function(x) rainbow_hcl(length(unique(x))))
    for(i in 1:length(newColsMicrobe)) {
      names(newColsMicrobe[[i]]) <- unique(microbe_meta[, i])
    }
    
    myColors <- append(newColsMeta, newColsMicrobe)
    
    return(
      pheatmap(
        microbe_x_metabolite,
        show_colnames = FALSE,
        show_rownames = FALSE,
        breaks = breaksList,
        treeheight_row = 0,
        treeheight_col = 0,
        fontsize_row = 6,
        color = diverge_hcl(length(breaksList)),
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        cluster_cols = TRUE,
        cluster_rows = TRUE,
        annotation_row = microbe_meta,
        annotation_col = metabolite_meta,
        # annotation_colors = myColors,
        cellwidth = 2,
        cellheight = 1,
        scale = "none"
      )
    )
  }

a <-
  prepare.blah(ihmp, microbe_abundance_table, metabolite_data, microbe_data)
b <- create.blah(a, microbe_data, metabolite_data)

ggsave(
  'microbeXmetabolite.png',
  plot = b,
  width = 7,
  height = 9,
  units = "in",
  dpi = 300
)
