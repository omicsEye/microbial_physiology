source('R/HeatMap.R')
library(RAM)

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


######################################################################
# TESTING ON iHMP DATA
# loading using load methods
microbe_data <- load.meta.data('Data/HMP/mxp_microbiome_v2019-06-25.csv')
microbe_data <- microbe_data[, c(7, 24, 25)] 
metabolite_data <- load.meta.data('Data/IHMP/HMDB_2019-07-12.csv', tax_column = 6)
metabolite_data <- metabolite_data[, c(91,92)]
sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv')
microbe_abundance_table <- load.abundance.data('Data/iHMP/iHMP_AbundanceTable_2019-07-09.csv')
microbe_abundance_table <- microbe_abundance_table[, -c(1)]

# read ihmp metabolite data
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_HILIC-neg_060517.csv',
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = TRUE,
  
  stringsAsFactors = FALSE
  
)

ihmp_col_names <- ihmp[4,-c(1:7)]

# filter out for only hmdb id
ihmp  <- ihmp[(ihmp$X.5) != "",]

# get column_names
ihmp_row_names <- ihmp[-1,6]

ihmp <- ihmp[-1,-c(1:7)]

ihmp <- apply(as.matrix(ihmp), 2, as.numeric)

# clean up ihmp data
colnames(ihmp) <- ihmp_col_names
rownames(ihmp) <- ihmp_row_names

# ajoined correlogram
a_table <- multi.correlogram(data_tables = list(ihmp, microbe_abundance_table), 
                             sample_datas = list(metabolite_data, microbe_data))

plot(a_table[[1]]$gtable)
plot(a_table[[2]]$gtable)
plot(a_table[[3]]$gtable)
plot(a_table[[4]]$gtable)

ihmp_map <- create.correlogram(ihmp, feature_meta = metabolite_data, omit = FALSE)
microbe_map <- create.correlogram(microbe_abundance_table, feature_meta = microbe_data, omit = FALSE)
