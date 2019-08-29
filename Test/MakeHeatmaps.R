source('R/HeatMap.R')
source('R/Utility.R')

# Microbe meta data
microbe_meta_data <- load.meta.data('Resources/Microbe Meta Tables/mxp_microbiome_v2019-08-26.csv')
microbe_meta_data <- microbe_meta_data[, c(7, 9, 16)]
microbe_meta_data <- clear.small.entry(microbe_meta_data, 0.05)

############################################################################
# CREATE HMP1-II MICROBE HEATMAPS
hmp_ot_abundance_table <- load.abundance.data('Resources/HMP/hmp1-II_metaphlan2-mtd-qcd.pcl_AbundanceTable_2019-07-09.csv')
hmp_ot_abundance_table <- hmp_ot_abundance_table[, -1]
hmp_ot_sample_data <- load.meta.data('Resources/HMP/SampleMetadata.csv')
hmp_ot_sample_data <- hmp_ot_sample_data[, c(1, 3)]

# heatmap for each morphological property
all.heatmap(data = hmp_ot_abundance_table,
            sample_meta = hmp_ot_sample_data,
            feature_meta = microbe_meta_data,
            filter = 'abundance')

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

save.figure(hmp_ot_correlogram, width = 35, height = 18)


############################################################################
# CREATE iHMP MICROBE HEATMAPS
ihmp_abundance_table <- load.abundance.data('Resources/iHMP/taxonomic_profiles.tsv_AbundanceTable_2019-07-09.csv')
ihmp_abundance_table <- ihmp_abundance_table[, -c(1)]
ihmp_microbe_sample_data <- load.meta.data('Resources/IHMP/hmp2_metadata.csv', tax_column = 2)
ihmp_microbe_sample_data <- ihmp_microbe_sample_data[, 70, drop = FALSE]

# heatmap for each morphological property
all.heatmap(data = ihmp_abundance_table,
            sample_meta = ihmp_microbe_sample_data,
            feature_meta = microbe_meta_data)

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

save.figure(ihmp_correlogram, width = 25, height = 10)

############################################################################
# CREATE iHMP METABOLITE HEATMAPS

# Metabolite meta data
metabolite_meta_data <- load.meta.data('Resources/Metabolite Meta Tables/HMDB_2019-07-30.csv')
metabolite_meta_data <- metabolite_meta_data[, c(11, 22, 23, 24, 25, 91:97), drop = FALSE]

ihmp_sample_data <- load.meta.data('Resources/IHMP/hmp2_metadata.csv', tax_column = 4)
ihmp_sample_data <- ihmp_sample_data[, 70, drop = FALSE]
##########
# C8-POS
ihmp_c8_intensity <- load.ihmp('Resources/iHMP/iHMP_metabolomics_C8-pos_060517.csv')
ihmp_c8_norm <- normalize.ihmp(ihmp_c8_intensity)

# heat map with all data
ihmp_c8_all_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = 0,
                                   show = FALSE, 
                                   omit_na = FALSE,
                                   unique_colors = FALSE)
save.figure(ihmp_c8_all_heatmap)
# heat map of 75th percentile based on abundance 
ihmp_c8_75a_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = .75,
                                   filter = 'abundance',
                                   show = FALSE, 
                                   omit_na = FALSE,
                                   unique_colors = FALSE)
save.figure(ihmp_c8_75a_heatmap)

# heat map of 75th percentile based on variance  
ihmp_c8_75v_heatmap <- create.heatmap(data = ihmp_c8_norm,
                                   sample_meta = ihmp_sample_data,
                                   feature_meta = metabolite_meta_data,
                                   percentile = .75,
                                   filter = 'variance',
                                   show = FALSE, 
                                   omit_na = FALSE,
                                   unique_colors = FALSE)
save.figure(ihmp_c8_75v_heatmap)

ihmp_c8_correlogram <- create.correlogram(
  data = ihmp_c8_intensity,
  feature_meta = metabolite_meta_data,
  show = TRUE,
  omit = FALSE,
  cluster_distance_method = 'euclidean',
  unique_colors = FALSE
)

save.figure(ihmp_c8_correlogram, width = 25, height = 15)

##########
# C18-NEG
ihmp_c18_intensity <- load.ihmp('Resources/iHMP/iHMP_metabolomics_C18-neg_060517.csv')
ihmp_c18_norm <- normalize.ihmp(ihmp_c18_intensity)

# heat map with all data
ihmp_c18_all_heatmap <- create.heatmap(data = ihmp_c18_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = 0,
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_c18_all_heatmap)
# heat map of 75th percentile based on abundance 
ihmp_c18_75a_heatmap <- create.heatmap(data = ihmp_c18_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = .75,
                                      filter = 'abundance',
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_c18_75a_heatmap)

# heat map of 75th percentile based on variance  
ihmp_c18_75v_heatmap <- create.heatmap(data = ihmp_c18_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = .75,
                                      filter = 'variance',
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_c18_75v_heatmap)

ihmp_c18_correlogram <- create.correlogram(
  data = ihmp_c18_intensity,
  feature_meta = metabolite_meta_data,
  show = TRUE,
  omit = FALSE,
  cluster_distance_method = 'euclidean',
  unique_colors = FALSE
)

save.figure(ihmp_c18_correlogram, width = 25, height = 15)

##########
# HILIC-POS
ihmp_hp_intensity <- load.ihmp('Resources/iHMP/iHMP_metabolomics_HILIC-pos_060517.csv')
ihmp_hp_norm <- normalize.ihmp(ihmp_hp_intensity)

# heat map with all data
ihmp_hp_all_heatmap <- create.heatmap(data = ihmp_hp_norm,
                                       sample_meta = ihmp_sample_data,
                                       feature_meta = metabolite_meta_data,
                                       percentile = 0,
                                       show = FALSE, 
                                       omit_na = FALSE,
                                       unique_colors = FALSE)
save.figure(ihmp_hp_all_heatmap, height = 12)
# heat map of 75th percentile based on abundance 
ihmp_hp_75a_heatmap <- create.heatmap(data = ihmp_hp_norm,
                                       sample_meta = ihmp_sample_data,
                                       feature_meta = metabolite_meta_data,
                                       percentile = .75,
                                       filter = 'abundance',
                                       show = FALSE, 
                                       omit_na = FALSE,
                                       unique_colors = FALSE)
save.figure(ihmp_hp_75a_heatmap, height = 12)

# heat map of 75th percentile based on variance  
ihmp_hp_75v_heatmap <- create.heatmap(data = ihmp_hp_norm,
                                       sample_meta = ihmp_sample_data,
                                       feature_meta = metabolite_meta_data,
                                       percentile = .75,
                                       filter = 'variance',
                                       show = FALSE, 
                                       omit_na = FALSE,
                                       unique_colors = FALSE)
save.figure(ihmp_hp_75v_heatmap, height = 12)

ihmp_hp_correlogram <- create.correlogram(
  data = ihmp_hp_intensity,
  feature_meta = metabolite_meta_data,
  show = TRUE,
  omit = FALSE,
  cluster_distance_method = 'euclidean',
  unique_colors = FALSE
)

save.figure(ihmp_hp_correlogram, width = 25, height = 15)

##########
# HILIC-NEG

ihmp_hn_intensity <- load.ihmp('Resources/iHMP/iHMP_metabolomics_HILIC-neg_060517.csv')
ihmp_hn_norm <- normalize.ihmp(ihmp_hn_intensity)

# heat map with all data
ihmp_hn_all_heatmap <- create.heatmap(data = ihmp_hn_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = 0,
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_hn_all_heatmap, height = 12)
# heat map of 75th percentile based on abundance 
ihmp_hn_75a_heatmap <- create.heatmap(data = ihmp_hn_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = .75,
                                      filter = 'abundance',
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_hn_75a_heatmap, height = 12)

# heat map of 75th percentile based on variance  
ihmp_hn_75v_heatmap <- create.heatmap(data = ihmp_hn_norm,
                                      sample_meta = ihmp_sample_data,
                                      feature_meta = metabolite_meta_data,
                                      percentile = .75,
                                      filter = 'variance',
                                      show = FALSE, 
                                      omit_na = FALSE,
                                      unique_colors = FALSE)
save.figure(ihmp_hn_75v_heatmap, height = 12)

ihmp_hn_correlogram <- create.correlogram(
  data = ihmp_hn_intensity,
  feature_meta = metabolite_meta_data,
  show = TRUE,
  omit = FALSE,
  cluster_distance_method = 'euclidean',
  unique_colors = FALSE
)

save.figure(ihmp_hn_correlogram, width = 25, height = 25)

############################################################################
# CREATE iHMP METABOLITE X MICROBE HEATMAPS

##########
# C8-POS
c8_m <- multi.correlogram(list(ihmp_abundance_table, ihmp_c8_norm), 
                  list(microbe_meta_data, metabolite_meta_data),
                  omit = FALSE)
save.figure(c8_m[[2]], width = 12, height = 20)
##########
# C18-NEG
c18_m <- multi.correlogram(list(ihmp_abundance_table, ihmp_c18_norm), 
                          list(microbe_meta_data, metabolite_meta_data))
save.figure(c18_m[[2]], width = 12, height = 20)
##########
# HILIC-POS
hp_m <- multi.correlogram(list(ihmp_abundance_table, ihmp_hp_norm), 
                          list(microbe_meta_data, metabolite_meta_data))
save.figure(hp_m[[2]], width = 12, height = 24)
##########
# HILIC-NEG
hn_m <- multi.correlogram(list(ihmp_abundance_table, ihmp_hn_norm), 
                          list(microbe_meta_data, metabolite_meta_data))
save.figure(hn_m[[2]], width = 12, height = 24)
