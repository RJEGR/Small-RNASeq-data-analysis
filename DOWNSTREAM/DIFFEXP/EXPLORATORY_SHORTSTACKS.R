# SHORSTACKS

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230308_test/"
list.files(path = path, pattern = "txt")

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)
res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

Count <- read_tsv(count_f)
Results <- read_tsv(res_f)

# Results %>% head() %>% view()
