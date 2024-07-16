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

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

f <- list.files(path = wd, pattern = "alignment_details.tsv", full.names = T)

alignment_details <- read_tsv(f)

readfile <- gsub(".clean.newid.subset.bam", "", basename(alignment_details$readfile))

alignment_details$readfile <- readfile

alignment_details <- alignment_details %>% mutate(rep = str_sub(readfile, end = -2))

alignment_details %>% group_by(readfile) %>% tally(count) # <- concordantly w/ initial proccessed reads

# 1)

ylab <- "Frac."

read_lengthL <- c("<21","21","22","23","24",">24")
repL <- c("HR248", "HR1108", "HR2476", "HR11076")

unique(alignment_details$rep)

# simplificar a mapeos unicos, multiples y sin mapeo

recode_to <- c(`U` = "Unique ", `P`= "Multiple",`H` = "Multiple", `R` = "Multiple", `N` = "Unmapped")

ylab <- "Reads (%)"
xlab <-  "Read Length"

getPalette <- c("#000056", "#2E71A7","#60A4CF", "#9ECAE1", "#EBE8E1") 

alignment_details %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  group_by(mapping_type, read_length) %>% 
  summarise(n = sum(count)) %>%
  group_by(read_length) %>% 
  pivot_wider(names_from = mapping_type, values_from = n ) %>% view()
  mutate(frac = n/sum(n)) %>% 
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) 
  ggplot(aes(x = read_length, y = frac, fill = mapping_type)) +
  geom_col(width = 0.8) + # position="dodge"
  # geom_segment(aes(xend = read_length, yend = 0,  color = mapping_type), linewidth = 1) +
  scale_y_continuous(ylab, labels = scales::percent) +
  labs(x = xlab) +
  # guides(fill = guide_legend(title = "Map")) + #  nrow = 5
  scale_fill_manual(values = c("#2E71A7", "#9ECAE1", "#EBE8E1")) +
  guides(fill = guide_legend(title = "", label_size = 0.5)) +
  # see::scale_fill_metro_d(reverse = T) +
  # see::scale_color_pizza(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 6.5) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    # legend.position = 'right',
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()) -> ps

ps <- ps + theme(  legend.position = 'top',
  legend.key.width = unit(0.15, "cm"),
  legend.key.height = unit(0.15, "cm"))

# https://easystats.github.io/see/articles/seecolorscales.html#overview-of-palette-colors

ggsave(ps, filename = 'ALIGNMENT_DETAILS_ES.png', path = path_out, width = 2, height = 1.25, device = png, dpi = 300)


# CALL FOR VARIANT AND PRECISION IDENTIFICATION:
which_vars <- c(c(as.character(21:24)))

Treads <- read_tsv(f) %>% select_at(all_of(which_vars)) %>% rowSums()

read_tsv(f) %>%
  mutate(PRECISION = Treads/Reads) %>%
  ggplot(aes(x = MIRNA, y = PRECISION)) +
  geom_boxplot()

which_vars <- c("Short", c(as.character(21:30)),  "Long")

sum(read_tsv(f)$Reads == Treads)

read_tsv(f) %>% 
  mutate(PRECISION = Treads/Reads) %>%
  select_at(all_of(c(which_vars, 'UniqueReads', 'MajorRNAReads', 'Reads', "FracTop", "PRECISION"))) %>%
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
  geom_histogram() + 
  facet_wrap(Strand ~ SRNAtype, scales = "free_y", nrow = 3, ncol = 3)

read_tsv(f) %>% 
  mutate(PRECISION = Treads/Reads) %>%
  ggplot(aes(PRECISION, FracTop)) + geom_point()

# FracTop: Fraction of Reads aligned to the top genomic strand.
# Strand: Inferred strandednes of the locus, inferred from FracTop and the --strand_cutoff setting