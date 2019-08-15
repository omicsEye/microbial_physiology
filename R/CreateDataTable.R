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

##########################################
# Source metabolite table
# HMDB
hmdb <- parse.hmdb()
