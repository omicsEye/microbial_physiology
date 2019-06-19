# m2interact
Extracts microbial physiological traits from online sources, merging them into a single, formated table

## Sourcing Microbial Data 

### Dependencies
Jupyter, RCurl, rjson, IRKernal

### Files
BacDiveApiCrawler.R
CleanProTrait.R
CombineData.R
SourceMicrobialData.ipynb

### Instructions
This script can be run using the provided Jupyter Notebook. Aside from the listed requirements, using BacDive requires an account associated with the site. This can be created at this [registration link](https://bacdive.dsmz.de/api/bacdive/registration/register/). After doing so, a username and password need to be provided to the script in order to access the API and make requests to it. 
Parameters can be supplied to each function, as well. Their options are listed with the function descriptions. 

### Functions
------------------------------------------------------------------------------------------------------------------------------
#### `BacDiveCrawler`

#### *Description*
BacDiveCrawler() retrieves information from the BacDive API, organizing it into a formatted table

#### *Usage*
`BacDiveCrawler(save_file = TRUE)`

#### *Arguments*
*save_file* if true, saves a .csv to the working directory containing the information extracted from the BacDive API

#### *Details*
Designed to traverse the API provided by BacDive. The BacDive API provides a database that can easily queried, providing microbial physiology data in the JSON format. Each specie contains its own ‘page’, which details information such as taxonomy, morphology, strain information, and more. This script currently selectively chooses certain traits to record, meaning that there is more data that could be chosen to extracted, if implemented. 

#### *Value*
Returns a data.frame containing information extracted from BacDive

#### *Warning*
Because this traverses the site’s API, it is still limited by internet speeds and the rate at which the site’s server responds. This can be detrimental to the speed at which the script can run. A possible future improvement would be to request multiple links at once, rather than one at a time. However, this, too, could prove problematic as the site’s server could be overloaded. 


#### `CleanProTrait`

#### *Description* 
CleanProTrait() retrieves information from a file downloaded from the ProTrait Atlas, formatting it into a table

#### *Usage* 
`CleanProTrait(save_file = TRUE)`

#### *Arguments* 
*save_file*	if true, saves a .csv to the working directory containing the information extracted from the ProTrait Atlas

#### *Details* 
Designed to extract information from a table created by ProTrait.  It lacks a format that generalizes traits, instead listing each type of trait (gram-positive, pathogenic in animals, aerobe, etc) as its column. Therefore, this script organizes this table into generalized traits, providing for an easy way to use this table for purposes such as annotation. 
It will first check if the ProTrait file already exists in the working directory. If it does not, it will download the file to the working directory and start formatting it. 

#### *Value* 
Returns a data.frame containing information extracted from ProTrait


#### `CombineData`

#### *Description*
CombineData() combines the given tables into a single, formatted table

#### *Usage*
`CombineData(protrait, bacdive, save_file = TRUE)` 

#### *Arguments*
*protrait* a data.frame containing the traits sourced from the ProTrait Atlas

*bacdive*	a data.frame containing the traits sources from the BacDive database

*save_file*	if true, saves a .csv to the working directory containing the information resulting from the combined table

#### *Details*
CombineData.R is a script that merges the tables extracted from BacDive and ProTrait. This script is required because the column labels produced for each of these tables are different and there are different traits extracted in general. This script works to create one cohesive table. 

#### *Value*
Returns a data.frame containing a information from the combined tables

#### *Warning*
This script lacks the functionality of being able to merge any two given tables. This therefore leaves it limited to the tables produced by the ProTrait and BacDive functions.
It also does not yet handle duplicate species. 

------------------------------------------------------------------------------------------------------------------------------

### Sources 
The [Bacterial Diversity Database](https://bacdive.dsmz.de), or “Bacdive”, is a database that provides microbial physiology data. It was created by and is maintained by the German Collection of Microorganisms and Cell Cultures GmbH. The [paper](https://academic.oup.com/nar/article/47/D1/D631/5106998) details how it was procured. 

The [ProTraits Atlas](http://protraits.irb.hr) describes microbial phenotypic traits. Multiple tables currently exist, and this script opts for data that is considered at least 90 % accurate in order to maximize the amount of accurate, recovered data. It was created by and is maintained by the Rudjer Boskovic Institute. Details about it can be found in the [paper](https://academic.oup.com/nar/article/44/21/10074/2290929). 
