source('R/HeatMap.R')

### FOR TESTING 
# load metadata
microbe_data <- read.delim(file.choose(), fill=NA, stringsAsFactors = FALSE, sep=',')
microbe_data <- microbe_data[!duplicated(microbe_data[,1]),]
row.names(microbe_data) <- microbe_data[,1]
species <- row.names(microbe_data)
microbe_data <- microbe_data[,-1]
# microbe_data <- microbe_data[, c(24, 25)]  # select certain traits
# microbe_data <- na.omit(microbe_data)


sample_data <- read.delim(file.choose(), sep = ",")
sample_data <- sample_data[!duplicated(sample_data[,1]),]
row.names(sample_data) <- sample_data[,1]
sample_data <- sample_data[,-c(1,3)]
# metabolite_data <- metabolite_data[,8:ncol(metabolite_data)]
# metabolite_data <- na.omit(metabolite_data)


abundance_table <- read.delim(file.choose(), sep = ",", row.names = 1)
species <- row.names(abundance_table)
abundance_table <- apply(abundance_table, 2, as.numeric)
row.names(abundance_table) <- species
samples <- colnames(abundance_table)
row.names(sample_data) <- samples

# data <- abundance_table
# row_meta <- microbe_data
# column_meta <- metabolite_data

# Create one v all heatmap
a <- one.v.all(abundance_table, sample_data, microbe_data, 
               percentile = 0.75, show = TRUE, 
               which = 2, column = 'Oxygen.tolerance', trait = 'anaerobe')

# Create heatmap
h <- create.heatmap(abundance_table, sample_data, microbe_data[, c(24, 7, 25)], show = TRUE, omit_na = FALSE)

# create all heatmaps
all.one.v.all(
  abundance_table, 
  sample_data, 
  microbe_data, 
  which=2, 
  column = 'met_util', 
  directory = 'Plots/')

# create correogram 
c <- create.correlogram(abundance_table, microbe_data[, c(24, 7, 25)], show = TRUE)

# create pcoa plot
ord_plots <- ordplots(abundance_table, as.data.frame(t(microbe_data[, 25, drop = FALSE])), output = '', outputname = 'pcoa' , method = 'pcoa')
