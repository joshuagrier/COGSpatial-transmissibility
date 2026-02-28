library(tidyverse)
library(rvest)

url <- "https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-vaccinations/covid-19-vaccinations-archive/covid-19-vaccinations-archive-weekly-publications/"
html <- read_html(url)


filevals <- html |> html_elements("p") |>
  html_elements("a") |>
  html_attr("href")

yy <- filevals[grepl("xlsx", filevals)]
xx <- grepl("weekly",yy)

library(HelpersMG)

for(f in yy[xx]) wget(f)
