# Access home page
# Access table
# Get all url on that page
# Get all data
# Move to next page
# Repeat

library(XML)
library(rvest)

url <- "http://bacmap.wishartlab.com/"
doc <- htmlParse(url)
links <- xpathSApply(doc, "//a/@href")

# free(doc)

which(grepl('/organisms/', links))

# check if can find next page link, if cannot or rejected, end program
which(grepl('/genomes?page=', links, fixed = TRUE))


url <- "http://bacmap.wishartlab.com/organisms/595"
meta_table <- url %>%
  read_html() %>%
  html_nodes(xpath='//*[@id="biography"]/table') %>%
  html_table(fill = TRUE)

url <- "http://bacmap.wishartlab.com/organisms/595"
a <- Async$new(urls = c(url))
b <- a$get()
c <- b[[1]]$parse("UTF-8")
d <- c %>%
  read_html() %>%
  html_nodes(xpath='//*[@id="biography"]/table') %>%
  html_table(fill = TRUE)

dtable <- meta_table[[1]]
dtable <- lapply(dtable, as.character)
dtable <- as.data.frame(dtable, stringsAsFactors = FALSE)
dtable <- t(dtable)
dtable <- as.data.frame(dtable, stringsAsFactors = FALSE)
names(dtable) <- dtable[1,]
dtable <- dtable[-1, , drop = FALSE]
