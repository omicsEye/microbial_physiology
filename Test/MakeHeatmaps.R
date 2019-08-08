source('R/HeatMap.R')
source('R/Utility.R')

# Microbe meta data
microbe_meta_data <- load.meta.data('mxp_microbiome_v2019-08-08.csv')
microbe_meta_data <- microbe_meta_data[, c(7, 9, 16)]
microbe_meta_data <- clear.small.entry(microbe_meta_data, 0.05)

############################################################################
# CREATE HMP1-II MICROBE HEATMAPS
hmp_ot_abundance_table <- load.abundance.data('Data/HMP/hmp1-II_metaphlan2-mtd-qcd.pcl_AbundanceTable_2019-07-09.csv')
hmp_ot_abundance_table <- hmp_ot_abundance_table[, -1]
hmp_ot_sample_data <- load.meta.data('Data/HMP/SampleMetadata.csv')
hmp_ot_sample_data <- hmp_ot_sample_data[, c(1, 3)]

# heat map with all data
hmp_ot_all_heatmap <- create.heatmap(data = hmp_ot_abundance_table,
                                 sample_meta = hmp_ot_sample_data,
                                 feature_meta = microbe_meta_data,
                                 percentile = 0,
                                 show = FALSE, 
                                 omit_na = FALSE)
save.figure(hmp_ot_all_heatmap)
# heat map of 75th percentile based on abundance 
hmp_ot_75a_heatmap <- create.heatmap(data = hmp_ot_abundance_table,
                                 sample_meta = hmp_ot_sample_data,
                                 feature_meta = microbe_meta_data,
                                 percentile = .75,
                                 filter = 'abundance',
                                 show = FALSE, 
                                 omit_na = FALSE)
save.figure(hmp_ot_75a_heatmap)

# heat map of 75th percentile based on variance  
hmp_ot_75v_heatmap <- create.heatmap(data = hmp_ot_abundance_table,
                                    sample_meta = hmp_ot_sample_data,
                                    feature_meta = microbe_meta_data,
                                    percentile = .75,
                                    filter = 'variance',
                                    show = FALSE, 
                                    omit_na = FALSE)
save.figure(hmp_ot_75v_heatmap)

hmp_ot_correlogram <- create.correlogram(
  data = hmp_ot_abundance_table,
  feature_meta = microbe_meta_data,
  show = FALSE,
  omit = FALSE,
  cluster_distance_method = 'euclidean'
)

ggsave(
  filename = paste0('', deparse(substitute(hmp_ot_correlogram)), '.png'),
  hmp_ot_correlogram,
  width = 35,
  height = 18,
  units = 'in',
  dpi = 300
)


############################################################################
# CREATE iHMP MICROBE HEATMAPS
ihmp_abundance_table <- load.abundance.data('Data/iHMP/taxonomic_profiles.tsv_AbundanceTable_2019-07-09.csv')
ihmp_abundance_table <- ihmp_abundance_table[, -c(1)]
ihmp_microbe_sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 2)
ihmp_microbe_sample_data <- ihmp_microbe_sample_data[, 70, drop = FALSE]

# heat map with all data
ihmp_all_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_microbe_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = 0,
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_all_heatmap)
# heat map of 75th percentile based on abundance 
ihmp_75a_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_microbe_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = .75,
                                     filter = 'abundance',
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_75a_heatmap)

# heat map of 75th percentile based on variance  
ihmp_75v_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_microbe_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = .75,
                                     filter = 'variance',
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_75v_heatmap)

ihmp_correlogram <- create.correlogram(
  data = ihmp_abundance_table,
  feature_meta = microbe_meta_data,
  show = FALSE,
  omit = FALSE,
  cluster_distance_method = 'euclidean'
)

ggsave(
  filename = paste0('', deparse(substitute(ihmp_correlogram)), '.png'),
  ihmp_correlogram,
  width = 25,
  height = 10,
  units = 'in',
  dpi = 300
)

############################################################################
# CREATE iHMP METABOLITE HEATMAPS

# Metabolite meta data
metabolite_meta_data <- load.meta.data('Data/iHMP/HMDB_2019-07-30.csv')
metabolite_meta_data <- metabolite_meta_data[, c(21, 22, 23, 25), drop = FALSE]

ihmp_sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 4)
ihmp_sample_data <- ihmp_sample_data[, 70, drop = FALSE]
##########
# C8-POS
ihmp_c8_intensity <- load.ihmp('Data/iHMP/iHMP_metabolomics_C8-pos_060517.csv')
ihmp_c8_norm <- normalize.ihmp(ihmp_c8_intensity)

# heat map with all data
ihmp_c8_all_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = 0,
                                   show = FALSE, 
                                   omit_na = FALSE)
save.figure(ihmp_c8_all_heatmap)
# heat map of 75th percentile based on abundance 
ihmp_c8_75a_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = .75,
                                   filter = 'abundance',
                                   show = FALSE, 
                                   omit_na = FALSE)
save.figure(ihmp_c8_75a_heatmap)

# heat map of 75th percentile based on variance  
ihmp_c8_75v_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = .75,
                                   filter = 'variance',
                                   show = FALSE, 
                                   omit_na = FALSE)
save.figure(ihmp_c8_75v_heatmap)

ihmp_c8_correlogram <- create.correlogram(
  data = ihmp_c8_intensity,
  feature_meta = metabolite_meta_data,
  show = FALSE,
  omit = FALSE,
  cluster_distance_method = 'euclidean'
)

ggsave(
  filename = paste0('', deparse(substitute(ihmp_c8_correlogram)), '.png'),
  ihmp_c8_correlogram,
  width = 25,
  height = 10,
  units = 'in',
  dpi = 300
)

##########
# C18-NEG
##########
# HILIC-POS
##########
# HILIC-NEG



############################################################################
# CREATE iHMP METABOLITE X MICROBE HEATMAPS

##########
# C8-POS
##########
# C18-NEG
##########
# HILIC-POS
##########
# HILIC-NEG
