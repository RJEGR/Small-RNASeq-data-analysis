
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

# USING SRNA INFO TO SPLIT:

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

DB %>% dplyr::count(SRNAtype)


DB <- DB %>% dplyr::select(Name, SRNAtype)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# 
pHpalette <- c(`Experimental`="#ad1f1f", `Control`= "#4575b4")
#

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

str(keep <- DB %>% filter(SRNAtype == "miR") %>% pull(Name))

sum(keep <- rownames(COUNTS) %in% keep)

nrow(COUNTS <- COUNTS[keep,])

# 2) BOXPLOT ====

recode_to <- c(`248` = "A) pH 8.0", `2476` = "B) pH 7.6", `1108`= "C) pH 8.0" ,`11076` = "D) pH 7.6")

recode_to2 <- c(`24 HPF` = "24 HPF", `110 HPF` = "110 HPF")

qprobs <- function(x) { 
  x <- x[x > 1]
  quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}

apply(log2(COUNTS+1), 2, qprobs) %>% 
  t() %>%
  as_tibble(rownames = 'LIBRARY_ID') %>%
  left_join(MTD) %>%
  mutate(sample_group = substr(LIBRARY_ID, 1,nchar(LIBRARY_ID)-1)) %>%
  dplyr::mutate(sample_group = dplyr::recode_factor(sample_group, !!!recode_to)) %>%
  mutate(sample_group = factor(sample_group, levels = recode_to)) %>%
  mutate(hpf = paste0(hpf, " HPF")) %>%
  mutate(hpf = factor(hpf, levels = recode_to2)) -> probs_df

probs_df %>%
  ggplot(., 
    aes(x = LIBRARY_ID, ymin = `5%`, lower = `25%`,
      middle = `50%`, upper = `75%`, ymax = `95%`)) +
  geom_errorbar(width = 0.3, position = position_dodge(0.6)) +
  geom_boxplot(aes(fill = pH), width = 0.5, stat = 'identity', 
    position = position_dodge(0.6)) +
  labs(y = expression('Expresión'~ '(' *log[2]* ')'), x = '') +
  theme_classic(base_family = "GillSans") +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme(
    legend.position = 'top',
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()) -> ptop

ptop <- ptop + 
  # facet_grid(~ sample_group, scales = 'free') +
  ggh4x::facet_nested(~ hpf+sample_group, nest_line = F, scales = "free") +
  theme(
    strip.background = element_rect(fill = 'grey89', color = 'white'))


# BOTTOM ====

apply(COUNTS, 2, function(x) sum(x > 0)) -> Total_genes

# Filter data by removing low-abundance genes

keep <- rowSums(edgeR::cpm(COUNTS) > 1) >= 2

nrow(count <- COUNTS[keep,])

# How singletones are per sample? ----

apply(count, 2, function(x) sum(x > 0)) -> filtered_genes

cbind(as_tibble(Total_genes, rownames = 'name'), as_tibble(filtered_genes)) -> n_genes

names(n_genes) <- c('LIBRARY_ID','Raw', 'Filt')

n_genes %>% mutate(pct_genes_retained = Filt/Raw) -> n_genes

n_genes 

n_genes <- n_genes %>%
  left_join(MTD) %>%
  mutate(sample_group = substr(LIBRARY_ID, 1,nchar(LIBRARY_ID)-1)) %>%
  dplyr::mutate(sample_group = dplyr::recode_factor(sample_group, !!!recode_to)) %>%
  mutate(sample_group = factor(sample_group, levels = recode_to)) 


n_genes %>%
  # mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = unique(LIBRARY_ID))) %>%
  ggplot() + 
  geom_col(aes(x = LIBRARY_ID, y = Raw, fill = pH)) + 
  # labs(x = 'Réplica', y = 'ARNs pequeños') +
  labs(x = 'Réplica', y = 'miRs') +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual("", values = rev(pHpalette)) +
  scale_fill_manual("", values = rev(pHpalette)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position = 'none',
    axis.text.x = element_text(angle = 45,
      hjust = 1, vjust = 1, size = 10)) -> pbottom

pbottom <- pbottom + 
  # facet_grid(~ sample_group, scales = 'free') +
  ggh4x::facet_nested(~ hpf+sample_group, nest_line = F, scales = "free") +
  theme(
    strip.background = element_blank(), 
    strip.text = element_blank())

library(patchwork)

ps <- ptop / plot_spacer() / pbottom + patchwork::plot_layout(heights = c(1,-0.25, 1.2))

ggsave(ps, filename = "transcripts_and_reads_plots_only_mirs.png", 
  path = path, width = 5, height = 4, device = png, dpi = 300)


view(n_genes)

view(data.frame(colSums(COUNTS)))
