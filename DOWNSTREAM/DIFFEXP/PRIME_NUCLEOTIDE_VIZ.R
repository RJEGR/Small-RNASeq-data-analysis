

# Q: VIZ 1ST NUCLEOTIE BIAS AT 5' ====
# THE 1ST NUC IS RELEVANT DUE TO THE AFFINITY OF AGOS TO CERTAIN CLASSES OF SRNAS ()
# Hoogstrate, S. W., ET AL (2014). https://doi.org/10.4161/worm.28234

# # The first base preference of known miRNA mature. 18~30-nt sRNAs were selected for analyzed and each histogram indicated the percentage of first base in the sRNAs with same RNA number ... Besides, we analyzed the first base preference for known mature miRNA (Fig. 2). Among these four groups: Ch-3, Ch-5, sample group A and sample group B showed preference towards U for the first base in 18~23-nt sRNAs, but the first base preference in 24~30-nt sRNAs was different between the four groups. sampleC and sampleD had more U preference than sampleA and samplegroupB in 24~35-nt sRNA

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

list.files(path = path, pattern = "txt")

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

RESULTS <- read_tsv(res_f)

SRNAS <- read_rds(paste0(path, "KNOWN_CLUSTERS_MIRS_PIRS.rds"))

str(which_pirs <- SRNAS %>% filter(grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())
str(which_mirs <- SRNAS %>% filter(!grepl("piR", KnownRNAs)) %>% distinct(Name) %>% pull())

str(query_names <- c(which_pirs, which_mirs))

# sum(!unique(RESULTS$Name) %in% query_names) # 65326

PRIME5_UNKNOWN <- RESULTS %>%
  filter(!Name %in% query_names) %>%
  group_by(MajorRNA, DicerCall) %>%
  summarise(MajorRNAReads = sum(MajorRNAReads), n = n()) %>%
  mutate(first_nuc = sapply(strsplit(MajorRNA, ""), `[`, 1)) %>%
  ungroup() %>%
  mutate(DicerCall = ifelse(DicerCall == "N", nchar(MajorRNA), DicerCall)) %>%
  mutate(biotype = "Unknown") %>%
  group_by(first_nuc, DicerCall, biotype) %>%
  summarise(MajorRNAReads = sum(MajorRNAReads), n = sum(n))

# PRIME5_UNKNOWN %>% view()

PRIME5_KNOWN <- RESULTS %>%
  # drop_na(KnownRNAs) %>% 
  filter(Name %in% query_names) %>%
  mutate(biotype = ifelse(Name %in% which_pirs, "piRs", "miRs")) %>%
  group_by(MajorRNA, biotype, DicerCall) %>%
  summarise(MajorRNAReads = sum(MajorRNAReads), n = n()) %>%
  mutate(first_nuc = sapply(strsplit(MajorRNA, ""), `[`, 1)) %>%
  ungroup() %>%
  mutate(DicerCall = ifelse(DicerCall == "N", nchar(MajorRNA), DicerCall)) %>%
  group_by(first_nuc, DicerCall, biotype) %>%
  summarise(MajorRNAReads = sum(MajorRNAReads), n = sum(n))


ylab <- "Number of reads"

rbind(PRIME5_KNOWN, PRIME5_UNKNOWN) %>%
  group_by(DicerCall, biotype) %>%
  # mutate(freq = MajorRNAReads/sum(MajorRNAReads)) %>%
  ggplot(aes(x = DicerCall, y = MajorRNAReads, fill = first_nuc)) + 
  facet_grid(biotype ~., scales = 'free') + 
  geom_col(width = 0.85) +
  see::scale_fill_social(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  scale_y_continuous(ylab, labels = scales::comma)
# scale_x_continuous(xlab, breaks = seq(min(out$Length),max(out$Length), by = 1)) +


# rbind(PRIME5_KNOWN, PRIME5_UNKNOWN) %>%
#   ggplot(aes(MajorRNAReads, n, color = biotype)) +
#   # facet_wrap(~ Length, scales = "free") +
#   geom_point(size = 5, alpha = 0.7) +
#   labs(caption = "Diff piR pattern expression from novel and miRs")


p1 <- rbind(PRIME5_KNOWN, PRIME5_UNKNOWN) %>%
  filter(biotype != "Unknown") %>%
  group_by(biotype) %>%
  mutate(MajorRNAReads = MajorRNAReads/sum(MajorRNAReads), n = n/sum(n)) %>%
  ggplot(aes(y = n, x = biotype, fill = first_nuc)) +
  geom_col(width = 0.85) +
  see::scale_fill_social(reverse = T) +
  labs(y = "Reads proprotion (%)") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  theme(legend.position = "top")



recode_to <- c(`MajorRNAReads` = "Number of reads ", `n`= "n clusters")

p2 <- rbind(PRIME5_KNOWN, PRIME5_UNKNOWN) %>%
  filter(biotype != "Unknown") %>%
  pivot_longer(cols = c("MajorRNAReads","n")) %>%
  dplyr::mutate(name = dplyr::recode_factor(name, !!!recode_to)) %>%
  ggplot(aes(y = value, x = biotype, fill = first_nuc)) +
  facet_grid(name ~., scales = "free_y") +
  geom_col(width = 0.85) +
  see::scale_fill_social(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 12) +
  scale_y_continuous("", labels = scales::comma) +
  theme(legend.position = "none")

p2

library(patchwork)

p1/p2
