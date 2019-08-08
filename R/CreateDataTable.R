source('R/CleanProTrait.R')
source('R/ParseIJSEM.R')
source('R/BacMapCrawler.R')
source('R/BacDiveApiCrawler.R')
source('R/ParseHMDB.R')
source('R/CombineData.R')

# Source microbe information
##########################################
# Protraits
protrait <- CleanProtrait()

# IJSEM
ijsem <- parse.ijsem()

# BacMap
bacmap <- bacmap.crawler()

# BacDive 
bacdive <- BacDiveCrawler()

##########################################
# Combine microbe tables

##########################################
# Source metabolite table
# HMDB
hmdb <- parse.hmdb()
