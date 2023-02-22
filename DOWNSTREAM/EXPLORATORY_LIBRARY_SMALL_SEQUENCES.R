
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(Biostrings)

library(tidyverse)

mtd <- read_tsv(list.files(path = '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/', pattern = 'METADATA', full.names = T))

path <- '~/Documents/MIRNA_HALIOTIS/RAW_DATA/usftp21.novogene.com/00.CleanData/'

subdir <- list.files(path = path, pattern = 'HR*', full.names = T)

fileName <- list.files(path = subdir, pattern = 'clean.fa.gz', full.names = T)

seqs1 <- readDNAStringSet(fileName[1], format="fasta")

count_widths <- function(fileName) {
  
  seqs <- readDNAStringSet(fileName, format="fasta")
  
  widths <- sort(width(seqs))
  widths <- as.data.frame(table(widths))
  widths$fileName <-  basename(fileName)
  
  widths <- as_tibble(widths)
  
  return(widths)
}

count_widths(fileName[1])

df <- lapply(fileName, count_widths)

data <- do.call(rbind, df)

data %>% mutate(fileName = sub('.clean.fa.gz', '', fileName)) %>% 
  left_join(mtd, by = c('fileName'='File')) %>% 
  ggplot(aes(x = widths, y = Freq)) + 
  # facet_grid(fileName ~ .) +
  facet_grid(hpf ~ .) +
  geom_col() +
  # geom_col(aes(fill = as.factor(Replicate)), position = position_dodge2()) + 
  theme_classic(base_family = "GillSans", base_size = 16) +
  labs(y = 'Freq (1E6)', x = 'Sequence Length', caption = '.clean.fa.gz Files')

# ShortRead::readFastq()

# READ FASTA (uncompressed) ----
# reads_collapsed.fa from mapper OR mirdeep2 module

path <- '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/TEST_09022023/'

fileName <- list.files(path = path, pattern = 'reads_collapsed.fa', full.names = T)

seqs1 <- readDNAStringSet(fileName, format="fasta")

seqs1 <- unique(sort(seqs1)) # deduplic

seqs1[sort(width(seqs1))]

widths <- width(seqs1)
widths <- as.data.frame(table(widths))

widths %>% ggplot(aes(x = widths, y = Freq/100000)) + geom_col() + 
  theme_classic(base_family = "GillSans", base_size = 16) +
  labs(y = 'Freq (1E6)', x = 'sRNA Length')
