# Get genesys information with Rvest
library(rvest)
url <- "https://www.genesys-pgr.org/data/view/383"
kwt1955 <- url %>%
  html() %>%
  html_nodes(xpath='//*[@id="content"]/div/div[2]/div[5]/table') %>%
               html_table()
             kwt1955 <- kwt1955[[1]]


url <- "https://www.genesys-pgr.org/data/view/384"
kwt1956 <- url %>%
  html() %>%
  html_nodes(xpath='//*[@id="content"]/div/div[2]/div[5]/table') %>%
  html_table()
kwt1956 <- kwt1956[[1]]


"https://www.genesys-pgr.org/data/view/236"