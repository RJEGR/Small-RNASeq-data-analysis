# PLOT alignment_details.tsv

# alignment details ----

# mapping_type
# U: Uniquely mapped (not a multimapper).
# P: Multimapper placed using the method set by option --mmap.
# R: Multimapper placed at random.
# H: Very highly multi-mapped read (>=50 hits).
# N: Unmapped reads.

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

recode_to <- c(`U` = "Uniquely ", `P`= "Multimapper (mmap)",`R` = "Random mmap", `H` = "Highly mmap (>=50 hits) ", `N` = "Unmapped")

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out//"

f <- list.files(path = wd, pattern = "alignment_details.tsv", full.names = T)

alignment_details <- read_tsv(f)

readfile <- gsub(".clean.newid.subset.bam", "", basename(alignment_details$readfile))

alignment_details$readfile <- readfile

alignment_details <- alignment_details %>% mutate(rep = str_sub(readfile, end = -2))

alignment_details %>% group_by(readfile) %>% tally(count) # <- concordantly w/ initial proccessed reads

# 1)

ylab <- "Frac"

read_lengthL <- c("<21","21","22","23","24",">24")
repL <- c("HR248", "HR1108", "HR2476", "HR11076")

unique(alignment_details$rep)

# 1) by replicate

alignment_details %>%
  group_by(mapping_type, read_length, readfile) %>% 
  summarise(n = sum(count)) %>%
  group_by(read_length, readfile) %>%
  mutate(frac = n/sum(n)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  mutate(rep = str_sub(readfile, end = -2)) %>%
  dplyr::mutate(rep = factor(rep, levels = repL)) %>%
  mutate(x_axis = str_sub(readfile, start = -1)) %>%
  ggplot(aes(x = x_axis, y = frac, fill = mapping_type)) +
  facet_grid(read_length ~ rep, scales = "free") +
  geom_col(width = 0.85) + # position="dodge"
  scale_y_continuous(ylab, labels = scales::percent) +
  labs(x = "Replicate") +
  guides(fill = guide_legend(title = "", nrow = 5)) +
  see::scale_fill_pizza(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'right',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0)) -> ps

ggsave(ps, filename = 'ALIGNMENT_DETAILS_FACET.png', path = wd, width = 6.7, height = 8, dpi = 300, device = png)

# by read size

# In concordance w/ spreadsheet sequencing

alignment_details %>%
  group_by(mapping_type,readfile) %>%
  summarise(n = sum(count)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  pivot_wider(names_from = mapping_type, values_from = n)

alignment_details %>%
  group_by(mapping_type,read_length) %>%
  summarise(n = sum(count)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  pivot_wider(names_from = mapping_type, values_from = n) 
  
# remark the mapping_type over k nucleotide read size

data <- alignment_details %>%
  group_by(mapping_type,read_length, readfile) %>%
  summarise(n = sum(count)) %>%
  mutate(ID = paste(mapping_type, read_length, sep = "-")) %>%
  ungroup() %>%
  select(-read_length, -mapping_type) %>%
  pivot_wider(names_from = readfile, values_from = n) 


data <- as(data[-1], "matrix")

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')


alignment_details %>%
  group_by(mapping_type,read_length, readfile) %>%
  summarise(n = sum(count)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  dplyr::mutate(mapping_type = factor(mapping_type, levels = recode_to)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  group_by(mapping_type, readfile) %>% 
  mutate(frac = n/sum(n)) %>% # tally(frac)
  ggplot(aes(y = readfile, x = read_length, fill = frac)) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_wrap(~ mapping_type, nrow = 1) +
  scale_fill_viridis_c(option = "B", name = "Frac.", direction = 1, 
    limits = c(0,1),
    labels = scales::percent_format(scale = 100)) +
  scale_x_discrete(position = 'top') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples) +
  theme_classic(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10)) -> p1

p1 <- p1 + theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  guides(fill = guide_colorbar(barheight = unit(0.2, "in"), 
    barwidth = unit(5, "in"), title.position = "left", label.position = "top",
    ticks.colour = "black", 
    frame.colour = "black",
    label.theme = element_text(size = 10, family = "GillSans")) )

# p1


p2 <- alignment_details %>%
  group_by(mapping_type,read_length) %>%
  summarise(n = sum(count)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  dplyr::mutate(mapping_type = factor(mapping_type, levels = recode_to)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  ggplot(aes(x = read_length, y = n)) +
  scale_x_discrete(position = 'top') +
  scale_y_reverse("Reads", labels = scales::comma) +
  facet_wrap(~mapping_type, nrow = 1, switch = "x") +
  geom_col(width = 0.85, fill = "black")  +
  labs(x = "") +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'right',
    panel.border = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()) 

library(patchwork)

ps <- p1/p2  + patchwork::plot_layout(heights = c(1,0.7))

ggsave(ps, filename = 'ALIGNMENT_DETAILS_HEATMAP.png', path = wd, width = 11, height = 6, dpi = 300, device = png)

# mapping_type read_length   count


alignment_details %>%
  group_by(mapping_type, read_length) %>% 
  summarise(n = sum(count)) %>%
  group_by(read_length) %>%
  mutate(frac = n/sum(n)) %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  ggplot(aes(x = read_length, y = frac, fill = mapping_type)) +
  geom_col(width = 0.85) + # position="dodge"
  scale_y_continuous(ylab, labels = scales::percent) +
  guides(fill = guide_legend(title = "")) + #  nrow = 5
  see::scale_fill_pizza(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'right',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) -> ps

# https://easystats.github.io/see/articles/seecolorscales.html#overview-of-palette-colors

ggsave(ps, filename = 'ALIGNMENT_DETAILS.png', path = wd, width = 5, height = 2.5, device = png)


# CHECK IF, DICERCALL PROFILE
# IT CORRESPOND TO MAJOR READ LENGTH THAT CONTRIBUTES TO CLUSTERING SNRAS
# If < 80% of the aligned reads are in the --dicermin and --dicermax boundaries, DicerCall is set to 'N'. Loci with a DicerCall of 'N' are unlikely to be small RNAs related to the Dicer-Like/Argonaute system of gene regulation.

f <- list.files(path = wd, pattern = "Results.txt", full.names = T)

m <- read_tsv(f) %>% 
  select(DicerCall, Short,c(as.character(21:30)),  Long)  

sample_cor = cor(as(m[-1], "matrix"), method='pearson', use='pairwise.complete.obs')

heatmap(sample_cor)

Levels_to <- c("<21",c(as.character(21:30)), ">30")

recode_to <- c(`Long` = ">30", `Short`= "<21")


m <- m %>% pivot_longer(-DicerCall) %>% 
  filter(value > 0) %>% 
  dplyr::mutate(name = dplyr::recode_factor(name, !!!recode_to)) %>%
  group_by(DicerCall,name) %>%
  summarise(value = sum(value)) 

x <- m %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  select(all_of(c("DicerCall",Levels_to)))

rowNames <- x$DicerCall

x <- as(x[-1], "matrix")

rownames(x) <- rowNames

sample_cor = cor(t(x), method='pearson', use='pairwise.complete.obs')
sample_dist = dist(x, method='euclidean')
hc_samples = hclust(sample_dist, method='complete')

heatmap(sample_cor)

sum(m$value) # 83,808,529

m %>%
  dplyr::mutate(name = factor(name, levels = Levels_to)) %>%
  dplyr::mutate(DicerCall = factor(DicerCall, levels = c("N",c(as.character(21:30))))) %>%
  ggplot(aes(y = name, x = DicerCall, fill = log10(value))) + 
  geom_tile(color = 'white', linewidth = 0.2) +
  # facet_wrap(~ mapping_type, nrow = 1) +
  scale_fill_viridis_c(option = "B", name = "x", direction = 1) +
  #   limits = c(0,1),
  #   labels = scales::percent_format(scale = 100)) +
  scale_x_discrete(position = 'top') +
  theme_classic(base_size = 10, base_family = "GillSans") +
  labs(y = 'Read Length') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10))


# CALL FOR VARIANT AND PRECISION IDENTIFICATION:
which_vars <- c(c(as.character(21:24)))

Treads <- read_tsv(f) %>% select_at(all_of(which_vars)) %>% rowSums()

read_tsv(f) %>%
  mutate(PRECISION = Treads/Reads) %>%
  ggplot(aes(x = MIRNA, y = PRECISION)) +
  geom_boxplot()

which_vars <- c("Short", c(as.character(21:30)),  "Long")

sum(read_tsv(f)$Reads == Treads)

read_tsv(f) %>% select_at(all_of(c(which_vars, 'UniqueReads', 'MajorRNAReads', 'Reads', "FracTop"))) %>%
  rstatix::cor_mat() %>%
  rstatix::cor_reorder() %>%
  rstatix::pull_lower_triangle() %>%
  rstatix::cor_plot(label = TRUE)

which_vars2 <- c('DicerCall','UniqueReads', 'MajorRNAReads', 'Reads')

# ???

read_tsv(f) %>% select_at(all_of(c(which_vars, which_vars2))) %>%
  mutate(PRECISION = (UniqueReads+MajorRNAReads)/Reads) %>%
  ggplot(aes(x = DicerCall, y = PRECISION)) +
  geom_boxplot()
  

DB %>% mutate(PRECISION = (UniqueReads+MajorRNAReads)/Reads) %>%
  ggplot(aes(fill = SRNAtype, PRECISION)) +
  geom_histogram() + facet_wrap(~ SRNAtype, scales = "free_y")

DB %>%
  ggplot(aes(fill = SRNAtype, FracTop)) +
  geom_histogram() + facet_wrap(~ SRNAtype, scales = "free_y")
