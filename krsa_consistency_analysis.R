# Explore the distribution and the correlation of the Z scores
# across experimental conditions.

suppressPackageStartupMessages({
    library(readr)
    library(purrr)
    library(dplyr)
})

creedenzymatic_files <- list.files("results")
