library(rBiopaxParser)

#' Draw enriched pathways
#' Visualize a metabolic pathway, highlighting metabolites that are enriched
#' @details plots the metabolic pathway
#' 
#' @param biopax_path the file path to a biopax representation of the metabolic pathway
#' @param metabolites a list of metabolites that are enriched
#' 
#' @return graphNEL representation of the metabolic pathway 
#' @export
draw.enriched.pathway <- function(biopax_path, metabolites) {
  biopax = readBiopax(biopax_path, verbose = FALSE)
  pw_list = listInstances(biopax, class="pathway")
  
  g <- NULL
  count <- 1
  while (is.null(g)) {
    g <- pathway2Graph(biopax, pw_list[[1]][count])
    count <- count + 1
  }
  # get list of node labels and find correct names
  node_names <- g@nodes
  p_datatable <-biopax$dt
  good_node_names <- c()
  
  for(i in node_names) {
    metabolite_local <- which(p_datatable$id == i)
    new_name <- paste(p_datatable$property_value[intersect(p_name_local, metabolite_local)], "\n", i)
    good_node_names <- append(good_node_names, new_name)
  }
  
  g@nodes <- good_node_names
  names(g@edgeL) <- good_node_names
  bio_nodes <- g@nodes
  
  nNodes <- list()
  
  enriched_nodes <- metabolites
  names(enriched_nodes) <- as.list(enriched_nodes)
  nNodes$label <- enriched_nodes
  
  color <-
    grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  label_colors <- color[seq(18, length(color), floor(length(color) / length(metabolites)))]
  
  names(label_colors) <- as.list(enriched_nodes)
  nNodes$color <- label_colors
  
  plot(g, nodeAttrs=nNodes, attrs=list(node=list(fontsize=60, cex=1)))  # neato, twopi
  return(g)
}
