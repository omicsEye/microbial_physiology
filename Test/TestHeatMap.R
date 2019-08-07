setwd("~/Documents/R_WorkSpace/m2interact")

source('R/HeatMap.R')
library(RAM)
library(recommenderlab)
source('R/Utility.R')

# TESTING ON iHMP DATA
# loading using load methods
microbe_data <- load.meta.data('Data/HMP/Final_v2019-07-30.csv')
microbe_data <- microbe_data[, c(7, 24, 25)] 
metabolite_data <- load.meta.data('Data/IHMP/HMDB_2019-07-30.csv', tax_column = 1)
metabolite_data <- metabolite_data[, c(11, 20:25, 91:97)]
sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 4)
sample_data <- sample_data[, 70, drop = FALSE]
microbe_sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 2)
microbe_sample_data <- microbe_sample_data[, 70, drop = FALSE]
microbe_abundance_table <- load.abundance.data('Data/iHMP/taxonomic_profiles.tsv_AbundanceTable_2019-07-09.csv')
microbe_abundance_table <- microbe_abundance_table[, -c(1)]

# read ihmp metabolite data
# iHMP_metabolomics_HILIC-neg_060517
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_HILIC-neg_060517.csv',
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = TRUE,
  
  stringsAsFactors = FALSE
  
)

# iHMP_metabolomics_HILIC-pos_060517
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_HILIC-pos_060517.csv',
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = TRUE,
  
  stringsAsFactors = FALSE
  
)

# read ihmp metabolite data
# iHMP_metabolomics_C8-pos_060517
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_C8-pos_060517.csv',
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = TRUE,
  
  stringsAsFactors = FALSE
  
)

# iHMP_metabolomics_C18-neg_060517
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_C18-neg_060517.csv',
  
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

ihmp[is.na(ihmp)] <- 0
p.1 <- t(apply(ihmp, 1, scale))
p.2 <- as(ihmp, "realRatingMatrix")
p.2 <- normalize(p.2, method="Z-score", row=TRUE)
p.2 <- as.matrix(p.2@data)
p.2[is.na(p.2)] <- 0

# ajoined correlogram
a_table <- multi.correlogram(data_tables = list(ihmp, microbe_abundance_table), 
                             sample_datas = list(metabolite_data, microbe_data))

plot(a_table[[1]]$gtable)
plot(a_table[[2]]$gtable)
plot(a_table[[3]]$gtable)
plot(a_table[[4]]$gtable)

ggsave(
  'microbeXmetaboliteExtra2.png',
  plot = a_table[[2]],
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

ihmp_map <-
  create.heatmap(
    p.2,
    sample_meta = sample_data,
    feature_meta = metabolite_data,
    omit = TRUE,
    percentile = 0,
    show = FALSE
  )

ggsave(
  'wow.png',
  plot = ihmp_map,
  width = 5,
  height = 6,
  units = "in",
  dpi = 300
)

d <- dist(p.2)
h <- hclust(d)


ihmp_norm_corr <- create.correlogram(p.2,
                   feature_meta = metabolite_data,
                   show = FALSE,
                   omit = FALSE)

ggsave(
  'normalized_correlation.png',
  plot = ihmp_norm_corr,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)

ihmp_nonnorm_corr<- create.correlogram(ihmp,
                   feature_meta = metabolite_data,
                   show = FALSE,
                   omit = FALSE)

ggsave(
  'nonnormalized_correlation.png',
  plot = ihmp_nonnorm_corr,
  width = 11,
  height = 7,
  units = "in",
  dpi = 300
)

all.heatmap(
  p.2,
  sample_meta = sample_data,
  feature_meta = metabolite_data,
  omit_na = TRUE,
  percentile = 0,
  show = FALSE
)

##########



m <- apply(microbe_data, c(1, 2), order.string)
m <- as.data.frame(m, stringsAsFactors = FALSE)

microbe_map <-
  create.heatmap(
    only_cd_abundance,
    sample_meta = microbe_sample_data,
    feature_meta = m,
    omit = FALSE,
    percentile = 0,
    show = TRUE,
    cluster_distance_method = "euclidean"
  )

ggsave(
  'microbe_heatmap.png',
  plot = microbe_map,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

microbe_corr <-
  create.correlogram(microbe_abundance_table,
                     feature_meta = microbe_data,
                     omit = FALSE,
                     cluster_distance_method = "euclidian")

ggsave(
  'microbe_correlogram.png',
  plot = microbe_corr,
  width = 30,
  height = 35,
  units = "in",
  dpi = 300
)

all.heatmap(
  microbe_abundance_table,
  sample_meta = microbe_sample_data,
  feature_meta = microbe_data,
  omit_na = TRUE,
  percentile = 0,
  show = FALSE
)

all.one.v.all(
  microbe_abundance_table,
  sample_meta = microbe_sample_data,
  feature_meta = m,
  column = 'diagnosis',
  percentile = 0.75,
  which = 1
)

