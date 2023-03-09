rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

mtd <- read_tsv('~/Documents/MIRNA_HALIOTIS/METADATA.tsv') %>% rename("sample_id" = "File")

path <- '~/Documents/MIRNA_HALIOTIS/PROFILING_BY_READ_LENGTH/'

f_list <- list.files(path = path, pattern = "clean.newid.profiling", full.names = T)

library(tidyverse)

read_df_and_summarise <- function(x) {
  
  recode_to <- c(`mirna` = "miRNAs", `unknown`= "Unknown",`rrna` = "rRNA", `trna` = "tRNA", `artifacts` = "Artifacts")
  
  sample_id <- sapply(strsplit(basename(x), "[.]"), `[`, 1)
  
  cat("\nReading sample: ", sample_id, "\n")
  
  # df <- read_tsv(x, col_names = T)
  
  .df <- data.table::fread(x)
  
  # seq_id <- sapply(strsplit(.df$`#name`, " "), `[`, 1)
  
  first_nuc <- sapply(strsplit(.df$seq, ""), `[`, 1)
  
  rnatype <- sapply(strsplit(.df$`#name`, " "), `[`, 2)
  
  clade <- sapply(strsplit(.df$`#name`, " "), `[`, 3)
  
  Length <- .df$length


  out <- data.frame(rnatype, Length, first_nuc) %>% count(rnatype, Length, first_nuc) 
  
  out <- dplyr::mutate(out, sample_id = sample_id)
  
  # Total <- nrow(.df)
  # out <- dplyr::mutate(out, pct = n/Total)
  
  out <- out %>%
    dplyr::as_tibble() %>% 
    dplyr::mutate(rnatype = gsub("rnatype:", "", rnatype)) %>%
    dplyr::mutate(rnatype = dplyr::recode_factor(rnatype, !!!recode_to))
  
  cat("\nSample", sample_id, "Done!!!\n")
  
  return(out)
  
}


# read_df_and_summarise(f_list[1])

# out <- lapply(f_list, read_df_and_summarise) # Take a while

# out <- do.call(rbind, out)

write_rds(out, file = paste0(path, "prof_by_read_length_summary.rds"))

out <- read_rds(paste0(path, "prof_by_read_length_summary.rds")) %>% left_join(mtd)

# Identification of known and novel miRNAs: 
#The first base preference of known miRNA mature. 18~30-nt sRNAs were selected for analyzed and each histogram indicated the percentage of first base in the sRNAs with same RNA number ... Besides, we analyzed the first base preference for known mature miRNA (Fig. 2). Among these four groups: Ch-3, Ch-5, sample group A and sample group B showed preference towards U for the first base in 18 ~ 23-nt sRNAs, but the first base preference in 24 ~ 30-nt sRNAs was different between the four groups. sampleC and sampleD had more U preference than sampleA and samplegroupB in 24 ~ 35-nt sRNA

xlab <- "Read Length (nt)"
ylab <- "Percentage of reads"

out %>% 
  mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  group_by(group, Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(group, Length) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  facet_grid(group ~ .) + geom_col() +
  scale_y_continuous(ylab, labels = scales::percent) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_pizza() +
  theme_bw(base_family = "GillSans") -> ps

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH.png', path = path, width = 12, height = 6.5)

# out %>% ggplot(aes(x = Length, y = n, fill = rnatype)) + geom_col()

# BY PCT

df %>% ggplot(aes(x = Length, y = pct, fill = first_nuc)) + geom_col() + 
  facet_wrap(~ rnatype, scales = 'free_y', nrow = 5)
