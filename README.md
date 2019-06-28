# m2interact
Extracts microbial physiological traits from online sources, merging them into a single, formated table

### Dependencies
Jupyter, RCurl, rjson, IRKernal

## Instructions

### Sourcing Microbial Data 
This script can be run using the provided Jupyter Notebook. Aside from the listed requirements, using BacDive requires an account associated with the site. This can be created at this [registration link](https://bacdive.dsmz.de/api/bacdive/registration/register/). After doing so, a username and password need to be provided to the script in order to access the API and make requests to it. 
Parameters can be supplied to each function, as well. Their options are listed with the function description found in the [wiki](https://github.com/broadinstitute/m2interact/wiki/Finding-Metadata-for-Microbial-Physiology-Traits/) 

### Creating Heatmaps
Abundance data should be loaded from a .csv. A path and the column that contains the feature labels must be supplied. 
`abundance_table <- load.abundance.data(path = <DATA>, column = <NUMBER>)`

Meta data should be loaded from a .csv. A path and the column that contains the feature labels must be supplied.
`meta_data <- load.meta.data(path = <DATA>, tax.column = <NUMBER>)`
This can be used to load both sample and feature metadata 

In order to create a heatmap, use the following function, supplying the abundance table, sample metadata, feature metadata. 
It includes the option of whether to show the heatmap when procuded, and whether to eliminate features that are missing metadata. 
`heatmap <- create.heatmap(abundance_table, sample_data, microbe_data[, c(24, 7, 25)], show = TRUE, omit_na = FALSE)`

------------------------------------------------------------------------------------------------------------------------------

### Sources 
The [Bacterial Diversity Database](https://bacdive.dsmz.de), or “Bacdive”, is a database that provides microbial physiology data. It was created by and is maintained by the German Collection of Microorganisms and Cell Cultures GmbH. The [paper](https://academic.oup.com/nar/article/47/D1/D631/5106998) details how it was procured. 

The [ProTraits Atlas](http://protraits.irb.hr) describes microbial phenotypic traits. Multiple tables currently exist, and this script opts for data that is considered at least 90 % accurate in order to maximize the amount of accurate, recovered data. It was created by and is maintained by the Rudjer Boskovic Institute. Details about it can be found in the [paper](https://academic.oup.com/nar/article/44/21/10074/2290929). 
