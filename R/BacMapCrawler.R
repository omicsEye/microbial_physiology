library(XML)
library(rvest)
library(crul)
library(data.table)
library(R.utils)

#' Get microbe urls
#' Retrive all microbial urls that exist on BacMap website
#'
#' @param url the url to a page containing a table of microbial urls, defaults to homepage
#' @param throttle the speed to ping the server
#'
#' @return list of urls links to microbial information
get.microbe.urls <-
  function(url = "http://bacmap.wishartlab.com/", throttle = 2) {
    microbe_urls <- list()
    next_page <- TRUE # whether there is a next page containing urls
    url <- url # current url on
    page <- 1 # which page expect to be on

    while (next_page) {
      start_time <- Sys.time() # for throttling

      doc <- htmlParse(url)

      # get all links that exist on page
      links <- xpathSApply(doc, "//a/@href")
      # get microbe urls
      microbe_urls <-
        append(microbe_urls, links[which(grepl('/organisms/', links))])

      # determine if next page
      next_url <- paste0('/genomes?page=', (page + 1))
      if (any(grepl(next_url, links, fixed = TRUE))) {
        url <- paste0("http://bacmap.wishartlab.com/", next_url)
        page <- page + 1
      } else {
        next_page <- FALSE
      }

      # can throttle in order to not ping server too often
      throttle_time <- throttle - (Sys.time() - start_time)
      if (throttle_time > 0)
        Sys.sleep(throttle_time)
    }

    return(unique(microbe_urls))
  }

#' Get microbe data
#' Get microbial phenotypic data from a url
#'
#' @param response the html response from the microbe url
#'
#' @return a dataframe of the mcirobial information from the html response
get.microbe.data <- function(response) {
  if (!is.atomic(response)) {
    response <- response$parse("UTF-8")
    if (is.character(response)) {
      meta_table <- response %>%
        read_html() %>%
        html_nodes(xpath = '//*[@id="biography"]/table') %>%
        html_table(fill = TRUE)

      dtable <- meta_table[[1]]
      dtable <- lapply(dtable, as.character)
      dtable <- as.data.frame(dtable, stringsAsFactors = FALSE)
      dtable <- t(dtable)
      dtable <- as.data.frame(dtable, stringsAsFactors = FALSE)
      names(dtable) <- dtable[1, ]
      dtable <- dtable[-1, , drop = FALSE]

      return(dtable)
    }

  }
  return(NA)
}

#' BacMap database crawler
#' Extract all microbial information from the BacMap website
#'
#' @param url the url of the bacmap homepage, defaults to known homepage
#' @param num_requests the number of requests to the website to perform at once
#'
#' @return a data frame containing all microbial information from the BacMap website
#' @export
bacmap.crawler <- function(url = "http://bacmap.wishartlab.com/", num_requests = 10) {
  microbe_urls <- get.microbe.urls()
  num_microbes <- length(microbe_urls)
  print(paste0('got ', num_microbes, ' urls'))
  current_page <- 1
  url <- url  # base url

  bacmap_table <- list()

  while (current_page <= num_microbes) {

    if (current_page > 50 && current_page %% 100 < 10) {
      num_results <- nrow(bacmap_table)
      message(paste("Gathered", num_results, "results"))  # prompt user about progress
    }

    # add urls to buffer
    url_buffer <- c()
    for (i in 1:num_requests) {
      url_buffer <-
        c(url_buffer,
          paste0(url, microbe_urls[[current_page]]))
      current_page <- current_page + 1

      if (current_page >= num_microbes) {
        break
      }
    }

    # get html
    responses <-
      try(withTimeout(Async$new(urls = url_buffer),
                      timeout = 5,
                      onTimeout = "error"))
    responses <- responses$get()

    # parse table
    if (class(responses) != "try-error") {
      for (i in 1:length(responses)) {
        current_row <- get.microbe.data(responses[[i]])

        if (is.data.frame(current_row)) {
          bacmap_table <-
            rbindlist(list(bacmap_table, current_row), fill = TRUE)
        }
      }
    }
  }

  bacmap_table <-
    apply(bacmap_table, 2, as.character)  # remove unwanted formating

  write.csv(bacmap_table,
            paste0("BacMap_v", Sys.Date(), ".csv"),
            row.names = FALSE)

  return(bacmap_table)
}
