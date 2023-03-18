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

# write_rds(out, file = paste0(path, "prof_by_read_length_summary.rds"))

out <- read_rds(paste0(path, "prof_by_read_length_summary.rds")) %>% left_join(mtd) %>%
  mutate(first_nuc = recode(first_nuc, `T` = "U"))

# 1) -----

xlab <- "Read Length (nt)"
ylab <- "Percentage of reads"

out %>% 
  mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  group_by(group, Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(group, Length) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  facet_grid(group ~ ., scales = "free_y") + geom_col(width = 0.85) +
  scale_y_continuous(ylab, labels = scales::percent) +
  # scale_y_continuous(ylab, labels = scales::comma) + # if y = n
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_pizza() -> bottom

bottom <- bottom + theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

# ggsave(bottom, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH.png', path = path, width = 6, height = 5, device = png)

# 1.1) top plot ---- 

xlab <- ""
ylab <- "Number of reads"

out %>% 
  mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  group_by(group, Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(group, Length) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  # facet_grid(group ~ ., scales = "free_y") + 
  geom_col(width = 0.85) +
  scale_y_continuous(ylab, labels = scales::comma) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_pizza() -> top

top <- top + theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank())

library(patchwork)

ps <- top / bottom + plot_layout(heights = c(0.3, 1), design = )

# t, b The top and bottom bounds of the area in the grid
# 
# l, r The left and right bounds of the area int the grid
# 

layout <- c(
  area(t = 0.5, b = 1,l = 1, r = 1),
  area(t = 1, b = 1, l = 1, r = 1)
)

top / bottom + plot_layout(heights = c(0.3, 1), design = layout)

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ.png', path = path, width = 6, height = 5, device = png)

# out %>% ggplot(aes(x = Length, y = n, fill = rnatype)) + geom_col()

# 2) BY PCT ----

# The first base preference of known miRNA mature. 18~30-nt sRNAs were selected for analyzed and each histogram indicated the percentage of first base in the sRNAs with same RNA number ... Besides, we analyzed the first base preference for known mature miRNA (Fig. 2). Among these four groups: Ch-3, Ch-5, sample group A and sample group B showed preference towards U for the first base in 18~23-nt sRNAs, but the first base preference in 24~30-nt sRNAs was different between the four groups. sampleC and sampleD had more U preference than sampleA and samplegroupB in 24~35-nt sRNA


xlab <- "Read Length (nt)"
ylab <- "Percentage of reads"

out %>% 
  # group_by(Length, rnatype, first_nuc) %>% summarise(n = sum(n)) %>%
  # group_by(Length, rnatype) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = n, fill = first_nuc)) + 
  facet_grid( rnatype ~., scales = 'free_y') + geom_col(width = 0.85) +
  # scale_y_continuous(ylab, labels = scales::percent) +
  scale_y_continuous(ylab, labels = scales::comma) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_social(reverse = T) -> ps

ps <- ps + theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH_NUC.png', path = path, width = 6, height = 5, device = png)


# 3) miRNA 1st position ----

ylab <- "Number of reads"
breaks_seq <- seq(18, 25, by = 1)

out %>% 
  filter(rnatype %in% 'miRNAs') %>%
  mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  group_by(Length, group, first_nuc) %>% summarise(n = sum(n)) %>%
  group_by(Length, group) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = n, fill = first_nuc)) + 
  facet_grid( group ~., scales = 'free_y') + geom_col(width = 0.85) +
  # scale_y_continuous(ylab, labels = scales::percent) +
  scale_y_continuous(ylab, labels = scales::comma) +
  scale_x_continuous(xlab, breaks = breaks_seq, limits = c(17, 26)) +
  see::scale_fill_social(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) -> ps

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH_NUC_kmiRNAs.png', path = path, width = 4, height = 4.5, device = png)
