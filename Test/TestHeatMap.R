setwd("~/Documents/R_WorkSpace/m2interact")

source('R/HeatMap.R')
library(RAM)

# TESTING ON iHMP DATA
# loading using load methods
microbe_data <- load.meta.data('Data/HMP/mxp_microbiome_v2019-06-25.csv')
# microbe_data <- microbe_data[, c(7, 24, 25)] 
metabolite_data <- load.meta.data('Data/IHMP/HMDB_2019-07-12.csv', tax_column = 6)
metabolite_data <- metabolite_data[, c(11, 12, 20:25, 92:95)]
sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 1)
microbe_sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 2)
microbe_sample_data <- microbe_sample_data[, c(75, 34, 40)]
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
    ihmp,
    sample_meta = sample_data,
    feature_meta = metabolite_data,
    omit = FALSE,
    percentile = 0
  )
ihmp_corr <-
  create.correlogram(
    ihmp,
    feature_meta = metabolite_data,
    sample_meta = sample_meta,
    show = TRUE,
    omit = FALSE
  )
all.heatmap(
  ihmp,
  sample_meta = sample_data,
  feature_meta = metabolite_data,
  omit_na = TRUE,
  percentile = 0,
  show = FALSE
)

##########

microbe_map <-
  create.heatmap(
    microbe_abundance_table,
    sample_meta = microbe_sample_data,
    feature_meta = microbe_data,
    omit = FALSE,
    percentile = 0,
    show = TRUE
  )

microbe_corr <-
  create.correlogram(microbe_abundance_table,
                     feature_meta = microbe_data,
                     omit = FALSE)
all.heatmap(
  microbe_abundance_table,
  sample_meta = microbe_sample_data,
  feature_meta = microbe_data,
  omit_na = TRUE,
  percentile = 0,
  show = FALSE
)

