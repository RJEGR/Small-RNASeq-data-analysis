rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

MTD <- read_tsv('~/Documents/MIRNA_HALIOTIS/METADATA.tsv')

scale_col <- c("#cd201f", "#FFFC00","#00b489","#31759b")

recode_to <- c(`Control` = "pH 8.0", `Low` = "pH 7.6")

MTD <- MTD %>%
  dplyr::rename("sample_id" = "LIBRARY_ID") %>%
  dplyr::mutate(hpf = paste0(hpf, " HPF")) %>%
  dplyr::mutate(hpf = factor(hpf, levels = c("24 HPF", "110 HPF"))) %>%
  dplyr::mutate(pH = dplyr::recode_factor(pH, !!!recode_to))

path <- '~/Documents/MIRNA_HALIOTIS/PROFILING_BY_READ_LENGTH/'

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

recode_to <- c(`miRNAs` = "miRNAs", `Unknown`= "Unknown",`rRNA` = "rRNA", `tRNA` = "tRNA", `Artifacts` = "Artifacts")

out <- read_rds(paste0(path, "prof_by_read_length_summary.rds")) %>% 
  left_join(MTD) %>%
  mutate(first_nuc = recode(first_nuc, `T` = "U")) %>%
  dplyr::mutate(rnatype = dplyr::recode_factor(rnatype, !!!recode_to))


# SANITY CHECK

out %>% group_by(sample_id) %>% tally(n)


# 1) Bottom -----

xlab <-  "Read Length (Nucleotides)" # "Read Length (nt)"
ylab <- "Sequences (Millions)"

recode_to <- c(`24 HPF`= "B) 24 HPF", `110 HPF` = "C) 110 HPF")

width_col <- 0.8
base_t_text_size <- 5

out %>% 
  dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  group_by(hpf, Length, rnatype) %>%
  summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>% 
  # summarise(n = sum(pct)) %>% # scheck
  ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  geom_col(width = width_col, color = "black", linewidth = 0.2) +
  # geom_segment(aes(xend = Length, yend = 0,  color = rnatype), linewidth = 1) +
  facet_grid(hpf ~ ., scales = "free_y") +
  scale_y_continuous(ylab, labels = scales::number_format(scale = 1/1000000, suffix = "")) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 4)) +
  see::scale_fill_pizza(reverse = T) +
  see::scale_color_pizza(reverse = T) -> bottom

bottom <- bottom + theme_classic(base_family = "GillSans", base_size = base_t_text_size) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.line.x = element_blank(),
    axis.text.x = element_text(size = base_t_text_size))

# 2) top -----

ylab <- "Frac. of reads"
xlab <- ""


top <- out %>% 
  group_by(Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>%
  mutate(sample_id = "A) Global") %>%
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  geom_col(width = width_col, color = "black", linewidth = 0.2) +
  # geom_segment(aes(xend = Length, yend = 0, color = rnatype), linewidth = 1) +
  facet_grid(sample_id ~ ., scales = "free_y") +
  scale_y_continuous(ylab, labels = scales::percent) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 2)) +
  see::scale_fill_pizza(reverse = T) + 
  see::scale_color_pizza(reverse = T) +
  guides(fill = guide_legend(title = "", label_size = 3)) +
  theme_classic(base_family = "GillSans", base_size = base_t_text_size) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank())


library(patchwork)


ps <- top / plot_spacer() / bottom + plot_layout(heights = c(0.5, -0.35, 1))

ggsave(ps, filename = 'FIGURE_2.png', path = path_out, width = 2, height = 1.7, device = png, dpi = 300)



# EXIT =======

out %>% 
  # mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  group_by(Length, rnatype) %>%
  summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>% 
  # summarise(n = sum(pct)) %>% # scheck
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  # facet_grid(group ~ ., scales = "free_y") + 
  geom_col(width = 0.85) +
  scale_y_continuous(ylab, labels = scales::percent) +
  # scale_y_continuous(ylab, labels = scales::comma) + # if y = n
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_pizza(reverse = T) -> top


top <- top + theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank())


library(patchwork)


ps <- top / plot_spacer() / bottom + plot_layout(heights = c(0.3, -0.1, 1))

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ.png', path = path, width = 6, height = 5.5, device = png)


# 1.1) top plot ---- 

xlab <- ""
# ylab <- "Número de lecturas"

out %>% 
  # mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  # group_by(group, Length, rnatype) %>% 
  # summarise(n = sum(n)) %>%
  # group_by(group, Length) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  # facet_grid(group ~ ., scales = "free_y") + 
  geom_col(width = 0.85) +
  scale_y_continuous(ylab, labels = scales::comma) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +
  see::scale_fill_pizza(reverse = T) -> top

top <- top + theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank())

library(patchwork)


ps <- top / plot_spacer() / bottom + plot_layout(heights = c(0.3, -0.1, 1))

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ.png', path = path, width = 6, height = 5.5, device = png)

# out %>% ggplot(aes(x = Length, y = n, fill = rnatype)) + geom_col()

# 2) BY PCT ----

# The first base preference of known miRNA mature. 18~30-nt sRNAs were selected for analyzed and each histogram indicated the percentage of first base in the sRNAs with same RNA number ... Besides, we analyzed the first base preference for known mature miRNA (Fig. 2). Among these four groups: Ch-3, Ch-5, sample group A and sample group B showed preference towards U for the first base in 18~23-nt sRNAs, but the first base preference in 24~30-nt sRNAs was different between the four groups. sampleC and sampleD had more U preference than sampleA and samplegroupB in 24~35-nt sRNA



xlab <-  "Read Length" # "Read Length (nt)"
ylab <- "Frecuency"

recode_to <- c(`miRs` = "A) miRs", `Desc`= "B) Desc.",`ARNr` = "C) ARNr",
  `ARNt` = "D) ARNt", `Artefactos` = "E) Artefactos")

table(out$rnatype)

  # 
  
out %>% 
  group_by(Length, rnatype, first_nuc) %>% summarise(n = sum(n)) %>%
  mutate(Freq = n / sum(n)) %>%
  dplyr::mutate(rnatype = dplyr::recode_factor(rnatype, !!!recode_to)) %>%
  # mutate(rnatype = factor(rnatype, levels = recode_to)) %>%
  ggplot(aes(x = Length, y = Freq, fill = first_nuc)) + 
  facet_grid( rnatype ~., scales = 'free_y') + geom_col(width = 0.85) +
  # scale_y_continuous(ylab, labels = scales::percent) +
  scale_y_continuous(ylab, labels = scales::comma) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 2)) +
  scale_fill_manual(values = scale_col) -> ps
  # see::scale_fill_bluebrown(reverse = T) -> ps

ps <- ps + 
  guides(fill = guide_legend(title = "")) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank())

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH_NUC.png', path = path, width = 4.5, height = 5, device = png, dpi = 300)


# 3) miRNA 1st position ----

ylab <- "Número de lecturas (Millones)"

breaks_seq <- seq(18, 26, by = 1)

out %>% 
  filter(grepl("miRs", rnatype)) %>%
  group_by(Length, hpf, pH, first_nuc) %>% summarise(n = sum(n)) %>%
  group_by(Length,  hpf, pH) %>% mutate(pct = n / sum(n)) %>%
  ggplot(aes(x = Length, y = n, fill = first_nuc)) + 
  # facet_grid( group ~., scales = 'free_y') + 
  ggh4x::facet_nested( hpf+pH ~ ., nest_line = F, scales = "free") +
  geom_col(width = 0.85) +
  # scale_y_continuous(ylab, labels = scales::percent) +
  # scale_y_continuous(ylab, labels = scales::comma) +
  scale_y_continuous(ylab, labels = scales::number_format(scale = 1/1000000, suffix = " M")) +
  scale_x_continuous(xlab, breaks = breaks_seq, limits = c(18, 26)) +
  # see::scale_fill_bluebrown(reverse = T) +
  scale_fill_manual("",values = scale_col) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) -> ps

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_LENGTH_NUC_kmiRNAs.png', path = path, width = 2.5, height = 3.5, device = png)

# BY:
recode_to <- c(`HR248` = "B)", `HR1108`= "C)",`HR2476` = "D)", `HR11076` = "E)")

# New facet label names for dose variable
y.labs <- c("B) pH 8.0", "C) pH 7.6", "D) pH 8.0", "E) pH 7.6")

names(y.labs) <- c("24 HPF:Ctrl pH", "24 HPF:Low pH", "110 HPF:Ctrl pH", "110 HPF:Low pH")

out %>% 
  group_by(hpf,pH, Length, rnatype) %>% 
  summarise(n = sum(n)) %>%
  mutate(pH = paste0(hpf, ":", pH)) %>%
  ggplot(aes(x = Length, y = n, fill = rnatype)) +
  ggh4x::facet_nested( hpf+pH ~ ., nest_line = F, scales = "free", labeller = labeller(pH = y.labs, .multi_line = T)) +
  # mutate(sample_id = substr(sample_id, 1,nchar(sample_id)-1)) %>%
  # mutate(group = paste0(hpf, "+", pH)) %>% # group = sample_id
  # group_by(sample_id, group, Length, rnatype) %>% summarise(n = sum(n)) %>%
  # dplyr::mutate(sample_id = dplyr::recode_factor(sample_id, !!!recode_to)) %>%
  # mutate(sample_id = factor(sample_id, levels = recode_to)) %>%
  # ggplot(aes(x = Length, y = n, fill = rnatype)) + 
  # facet_grid(sample_id ~ ., scales = "free_y") +
  geom_col(width = 0.85)  +
  # scale_y_continuous("Número de lecturas", labels = scales::comma) +
  scale_y_continuous("Número de lecturas", labels = scales::number_format(scale = 1/1000000, suffix = " M")) +
  scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 2)) +
  see::scale_fill_pizza(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'none',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) -> bottom

top <- out %>% 
  group_by(Length, rnatype) %>% summarise(n = sum(n)) %>%
  group_by(Length) %>% mutate(pct = n / sum(n)) %>%
  mutate(sample_id = "A) Global") %>%
  ggplot(aes(x = Length, y = pct, fill = rnatype)) + 
  geom_col(width = 0.85) +
  facet_grid(sample_id ~ ., scales = "free_y") +
  scale_y_continuous("Frecuencia", labels = scales::percent) +
  scale_x_continuous("", breaks = seq(min(out$Length),max(out$Length), by = 2)) +
  see::scale_fill_pizza(reverse = T) + 
  guides(fill = guide_legend(title = "")) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank())

ps <- top / bottom + plot_layout(heights = c(0.3, 1))

ps <- top / plot_spacer() / bottom + plot_layout(heights = c(0.3, -0.1, 1))

ggsave(ps, filename = 'PROFILING_SAMPLES_BY_READ_2.png', path = path, width = 4, height = 4.5, device = png, dpi = 300)

# PLOT CLADE SPECIFIC MIRS? ====

path <- '~/Documents/MIRNA_HALIOTIS/RAW_INPUT/MIRTRACE/mirtrace.20230107-184757.848/'

f <- list.files(path = path, pattern = "mirtrace-stats-contaminatio", full.names = T)

df <- read_tsv(f[1])

which_names <- df %>% select(starts_with("HR")) %>% names()

recode_to <- c(`HR248` = "A)", `HR1108`= "B)",`HR2476` = "C)", `HR11076` = "D)")

out <- df %>% 
  pivot_longer(cols = all_of(which_names), names_to = "sample_id", values_to = "Reads") %>% 
  left_join(MTD) %>%
  mutate(sample_group = substr(sample_id, 1,nchar(sample_id)-1)) %>%
  dplyr::mutate(sample_group = dplyr::recode_factor(sample_group, !!!recode_to)) %>%
  mutate(sample_group = factor(sample_group, levels = recode_to)) %>%
  filter(Reads > 0) %>%
  group_by(sample_id) %>%
  mutate(pct = Reads / sum(Reads)) %>%
  mutate(CLADE = stringr::str_to_sentence(CLADE))

CLADE_L <- out %>% group_by(CLADE) %>% tally(Reads, sort = T) %>% pull(CLADE)

ps <- out %>%
  mutate(sample_id = gsub("^HR", "", sample_id)) %>%
  mutate(CLADE = factor(CLADE, levels = CLADE_L)) %>%
  ggplot(aes(y = pct, x = sample_id, fill = CLADE)) + 
  facet_grid( ~ sample_group, scales = "free", space = "free") +
  geom_col() +
  labs(x = "", caption = "miR clado-especifico") +
  scale_y_continuous("% Lecturas", labels = scales::percent) +
  see::scale_fill_bluebrown(reverse = F) + 
  guides(fill = guide_legend(title = "", nrow = 1)) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank())
    


ggsave(ps, filename = 'CLADE_SPECIFIC_MIR.png', path = path, width = 7.2, height = 4, device = png)


out %>% group_by(CLADE, sample_id) %>% tally(pct, sort = T) %>%
  pivot_wider(names_from = sample_id, values_from = n)

out %>% 
  mutate(CLADE = factor(CLADE, levels = CLADE_L)) %>%
  mutate(sample_id = gsub("^HR", "", sample_id)) %>%
  ggplot(aes(x = sample_id, y = Reads)) + 
  facet_grid( ~ sample_group, scales = "free", space = "free") +
  geom_col() +
  scale_y_continuous("", labels = scales::comma)


out %>%
  mutate(CLADE = factor(CLADE, levels = CLADE_L)) %>%
  ggplot(aes(y = CLADE, x = sample_id, fill = pct)) + 
  # facet_grid( ~ sample_group, scales = "free") +
  geom_tile() +
  # scale_fill_grey()
  scale_fill_viridis_c(option = "B", direction = 1, 
    limits = c(0,1),
    labels = scales::percent_format(scale = 100)) +
  guides(fill = guide_legend(title = "", nrow = 1)) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.major = element_blank())

# read_tsv(f[2])
# 

# 
read_tsv(f[2]) %>% pivot_longer(cols = all_of(which_names)) %>% group_by(CLADE) %>%
  summarise(n = n(), Reads = sum(value))
