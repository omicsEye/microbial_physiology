# m2interact
**m2interact** is a physiology trait collection and data representation tool. It can be used to parse microbial and metabolite physiology data from online sources and construct a single data table for each. This datatable, combined with sample omics data (ie. microbial abundance tables & metabolite intensity values) can be represented using different visualization methods. 

## Instructions  
For a more indepth look at each function, refer to the [wiki](https://github.com/broadinstitute/m2interact/wiki)

### Sourcing Microbial Data 
Install the required dependencies using `install.packages(c(Jupyter, RCurl, rjson, IRKernal, pheatmap, ggplot2, RColorBrewer, XML, foreach, parallel, doParallel, data.table, utils, rlist, crul, jsonlite, R.utils, rvest, colorspace, recommenderlab, RAM))`

Be sure to set the working directory to the m2Interact folder using `setwd("~/<PATH>/m2Interact")`\
Source all files found in the *R* folder (if the package is not built) 

Using the BacDive API requires an account associated with the site. This can be created at this [registration link](https://bacdive.dsmz.de/api/bacdive/registration/register/). After doing so, a username and password need to be provided to the script in order to access the API and make requests to it. 
Parameters can be supplied to each function, as well.  

In order to create the table using default settings, run:
```R

# Source microbe information
##########################################
# Protraits
protrait <- clean.protrait()

# IJSEM
ijsem <- parse.ijsem()

# BacMap
bacmap <- bacmap.crawler()

# BacDive 
bacdive <- bacdive.crawler(<USERNAME>, <PASSWORD>)

##########################################
# Combine microbe tables
all_table <- list(protrait, ijsem, bacmap, bacdive)
total_table <- combine.data(all_table)
```
A .csv will be saved locally 

### Sourcing Metabolite Data
To parse the hmdb database, run:
```R
parse.hmdb()
```
A .csv will be saved locally

### Creating Heatmaps
An example of creating a heat map of hmp1-II data annotating it with microbioal physiology
```R
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

# heat map with all data
hmp_ot_all_heatmap <- create.heatmap(data = hmp_ot_abundance_table,
                                 sample_meta = hmp_ot_sample_data,
                                 feature_meta = microbe_meta_data,
                                 percentile = 0,
                                 show = FALSE, 
                                 omit_na = FALSE)
save.figure(hmp_ot_all_heatmap)
```



------------------------------------------------------------------------------------------------------------------------------

### Sources 
The [Bacterial Diversity Database](https://bacdive.dsmz.de), or “Bacdive”, is a database that provides bacterial physiology data. It was created by and is maintained by the German Collection of Microorganisms and Cell Cultures GmbH. The [paper](https://academic.oup.com/nar/article/47/D1/D631/5106998) details how it was procured. 

The [ProTraits Atlas](http://protraits.irb.hr) describes microbial phenotypic traits. Multiple tables currently exist, and this script opts for data that is considered at least 90 % accurate in order to maximize the amount of accurate, recovered data. It was created by and is maintained by the Rudjer Boskovic Institute. Details about it can be found in the [paper](https://academic.oup.com/nar/article/44/21/10074/2290929). 

The [BacMap Database](http://bacmap.wishartlab.com) provides bacterial phenotypic traits. It is developed and maintained by the Wishart Lab at the University of Alberta. More details are found in its [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245156/). 

The [database](https://figshare.com/articles/International_Journal_of_Systematic_and_Evolutionary_Microbiology_IJSEM_phenotypic_database/4272392) compiled by Albert Barberan from the [International Journal of Systematic and Evolutionary Microbiology](https://microbiologysociety.org) provides microbial phenotypic data. Details can be found in this [paper](https://msphere.asm.org/content/2/4/e00237-17). 

The [Human Metabolome Database](http://www.hmdb.ca) gives a combination of clinical, chemical, and biological properties of metabolites found in the human body. It was created and is maintained by the Wishart Lab as part of the Metabolomics Inovation Center. This [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753273/) gives more details. 
