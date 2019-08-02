source('R/HeatMap.R')
source('R/Utility.R')
source('Test/TestCleanTable.R')

# data tables
###
# Microbe meta data
microbe_meta_data <- load.meta.data('Data/HMP/Final_v2019-07-30.csv')
microbe_meta_data <- microbe_meta_data[, c(7, 24, 25)] 
microbe_meta_data <- apply(microbe_meta_data, c(1, 2), order.string)
microbe_meta_data <- as.data.frame(microbe_meta_data, stringsAsFactors = FALSE)
microbe_meta_data <- clear.small.entry(microbe_meta_data, 0.05)

# Metabolite meta data
metabolite_meta_data <- load.meta.data('Data/iHMP/HMDB_2019-07-30.csv')


############################################################################
# CREATE HMP1-II MICROBE HEATMAPS
hmp_ot_abundance_table <- load.abundance.data('Data/HMP/AbundanceTable.csv')
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



############################################################################
# CREATE iHMP MICROBE HEATMAPS
ihmp_abundance_table <- load.abundance.data('Data/iHMP/taxonomic_profiles.tsv_AbundanceTable_2019-07-09.csv')
ihmp_abundance_table <- ihmp_abundance_table[, -c(1)]
ihmp_sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 2)
ihmp_sample_data <- ihmp_sample_data[, 70, drop = FALSE]

# heat map with all data
ihmp_all_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = 0,
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_all_heatmap)
# heat map of 75th percentile based on abundance 
ihmp_75a_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = .75,
                                     filter = 'abundance',
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_75a_heatmap)

# heat map of 75th percentile based on variance  
ihmp_75v_heatmap <- create.heatmap(data = ihmp_abundance_table,
                                     sample_meta = ihmp_sample_data,
                                     feature_meta = microbe_meta_data,
                                     percentile = .75,
                                     filter = 'variance',
                                     show = FALSE, 
                                     omit_na = FALSE)
save.figure(ihmp_75v_heatmap)

############################################################################
# CREATE iHMP METABOLITE HEATMAPS

##########
# C8-POS
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
