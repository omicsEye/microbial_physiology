source('R/HeatMap.R')
library(massPattern)

# read ihmp data
ihmp <- read.csv(
  
  file.choose(),
  
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


p <- na.omit(ihmp)
total_col <- apply(p, 2, sum)

i <- 1
blah <- function(x) {
  a <- x / total_col[i]
  i <<- i + 1
  return(a)
}

p <- apply(p, 2, blah)

# get metabolite metadata
meta_data <- read.csv(
  
  file.choose(),
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE,
)
meta_data <- meta_data[,-c(1)]
meta_data <- meta_data[!duplicated(meta_data[,2]),]
meta_data <- meta_data[!is.na(meta_data$OldHMDBIDs), ]

row.names(meta_data) <- meta_data[,2]
meta_data <- meta_data[,-c(1:8)]

sample_data <- read.csv(
  
  file.choose(),
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE,
)
sample_data <- subset(sample_data,data_type == 'metabolomics')
rownames(sample_data) <- sample_data$site_sub_coll
sample_data <- sample_data[-1,-c(1:6)]
sample_data <- sample_data[, 'site_name', drop = FALSE]

create.heatmap(data = ihmp, sample_meta = sample_data, feature_meta = meta_data, show = TRUE)

create.correlogram(p, meta_data)

####################################################################

# metadata <- read.csv(
#   
#   file.choose(),
#   
#   header = TRUE,
#   
#   fill = FALSE,
#   
#   comment.char = "" ,
#   
#   check.names = TRUE
#   
# )
# 
data <- massPattern::load_data(

  input = file.choose(),

  type = 'all'

)

abundance_data <- data$data
abundance_data <- abundance_data[order(rownames(abundance_data)),]

# sample info
sample_info <- data$sample_metadata

# feature info ( e.g. m/z and RT)
features_info <- data$feature_metadata


ord_plots <- ordplots(abundance_data, sample_info, output = 'Plots/', outputname ='iHMP_PCOA_Feature_Test', method = 'pcoa')

# 
# metadata <- subset(metadata, data_type == 'metabolomics')
# rownames(metadata) <- metadata$Tube.A..Metabolomics
# metadata <- metadata[order(rownames(metadata)),]
# 
# metabolite_data <- read.csv(
#   
#   file.choose(),
#   
#   header = TRUE,
#   
#   fill = FALSE,
#   
#   comment.char = "" ,
#   
#   check.names = TRUE
#   
# ) 
# 
# rownames(metabolite_data) <- metabolite_data$Name
# metabolite_data <- metabolite_data[order(rownames(metabolite_data)),]
