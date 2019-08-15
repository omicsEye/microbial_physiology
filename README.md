# m2interact
**m2interact** is a physiology trait collection and data representation tool. It can be used to parse microbial and metabolite physiology data from online sources and construct a single data table for each. This datatable, combined with sample omics data (ie. microbial abundance tables & metabolite intensity values) can be represented using different visualization methods. 

## Instructions  
For a more indepth look at each function, refer to the [wiki](https://github.com/broadinstitute/m2interact/wiki)

### Sourcing Microbial Data 
Aside from the listed requirements, using BacDive requires an account associated with the site. This can be created at this [registration link](https://bacdive.dsmz.de/api/bacdive/registration/register/). After doing so, a username and password need to be provided to the script in order to access the API and make requests to it. 
Parameters can be supplied to each function, as well.  

To run it manually, execute: (This is also found in the R/CreateDataTable.R file)
```R
source('R/CleanProTrait.R')
source('R/ParseIJSEM.R')
source('R/BacMapCrawler.R')
source('R/BacDiveApiCrawler.R')
source('R/ParseHMDB.R')
source('R/CombineData.R')

# Source microbe information
##########################################
# Protraits
protrait <- clean.protrait()

# IJSEM
ijsem <- parse.ijsem()

# BacMap
bacmap <- bacmap.crawler()

# BacDive 
bacdive <- bacdive.crawler()

##########################################
# Combine microbe tables
all_table <- list(protrait, ijsem, bacmap, bacdive)
total_table <- combine.data(all_table)
```

### Sourcing Metabolite Data
To parse the hmdb database, run:
```R
parse.hmdb()
```
A file can be supplied which represents the name of the metabolite datatable file to parse (does not need to be provided if a local copy of the hmdb metabolite file does not exist). By default, it is *hmdb_metabolites.xml*.
A link can be supplied which represents the URL download link for the metabolite datatable. By default, it is *http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip*

### Creating Heatmaps
First, run, `source('R/HeatMap.R')`

Abundance data should be loaded from a .csv. A path and the column that contains the feature labels must be supplied.
NOTE: The path must be from the working directory
```R
abundance_table <- load.abundance.data(path = <DATA>, column = <NUMBER>)
```

Meta data should be loaded from a .csv. A path and the column that contains the feature labels must be supplied.
NOTE: The path must be from the working directory
```R
meta_data <- load.meta.data(path = <DATA>, tax.column = <NUMBER>)
```
This can be used to load both sample and feature metadata 

In order to create a heatmap, use the following function, supplying the abundance table, sample metadata, feature metadata. 
It includes the optional parameters of whether to filter clusters for a certain percentile, to show the heatmap when procuded, and whether to eliminate features that are missing metadata. 
```R
heatmap <- create.heatmap(data = <ABUNDANCE>, sample_meta = <SAMPLE>, feature_meta = <FEATURE>)
```

There is also the choice of creating heatmaps that compare only one feature type against all others in a feature category (ex. aerobic respiration v all other oxygen requirements). This requires the following function, with the new parameters being 
* `which` represents whether to use the sample(1) or feature(2) metadata
* `column` represents the name of the column of the feature category (ex. 'gram type')
* `trait` represents the feature trait to isoloate (ex. 'variable') 
```R
one_heatmap <- one.v.all(data = <ABUNDANCE>, sample_meta = <SAMPLE>, feature_meta = <FEATURE>, which = <NUMBER>, column = <CATEGORY>, trait = <TYPE>)
```

The following method can be used to repeat this 'one v all' process for every trait type in a trait category, creating a heatmap for each. The produced heatmaps can be saved to a desired directory. 
```R
all.one.v.all(data = <ABUNDANCE>, sample_meta = <SAMPLE>, feature_meta = <FEATURE>, which = <NUMBER>, column = <CATEGORY>, directory = <PATH>)
```

Instead of a sample x feature heatmap, a correlogram heatmap can be produced. Providing metadata and an abundance data, use the following function
```R
correlogram <- create.correlogram(data = <ABUNDANCE>, feature_meta = <METADATA>, show = TRUE)
```

------------------------------------------------------------------------------------------------------------------------------

### Sources 
The [Bacterial Diversity Database](https://bacdive.dsmz.de), or “Bacdive”, is a database that provides microbial physiology data. It was created by and is maintained by the German Collection of Microorganisms and Cell Cultures GmbH. The [paper](https://academic.oup.com/nar/article/47/D1/D631/5106998) details how it was procured. 

The [ProTraits Atlas](http://protraits.irb.hr) describes microbial phenotypic traits. Multiple tables currently exist, and this script opts for data that is considered at least 90 % accurate in order to maximize the amount of accurate, recovered data. It was created by and is maintained by the Rudjer Boskovic Institute. Details about it can be found in the [paper](https://academic.oup.com/nar/article/44/21/10074/2290929). 

The [Human Metabolome Database](http://www.hmdb.ca) gives a combination of clinical, chemical, and biological properties of metabolites found in the human body. It was created and is maintained by the Wishart Lab as part of the Metabolomics Inovation Center. This [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753273/) gives more details. 
