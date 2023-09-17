# TARGETSCAN
# ESTE ANALISIS ES DIRIGIDO A SOLO A LOS 131 MIRS:

library(tidyverse)


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

f <- "mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out"

f <- list.files(path = wd, pattern = f, full.names = T)

# Flat GTF info:

# utr_f <- list.files(path = wd, full.names = T, pattern = "three_prime_utr.ids")

# str(x <- read_lines(utr_f)) # 54 432 <-- i.e the N sequences from three_prime_utr.fa

# str(target <- sapply(strsplit(x, " "), `[`, 1))

# str(gene_id <- sapply(strsplit(x, " "), `[`, 2)) 

# nrow(utr_source <- data.frame(target, gene_id) %>% as_tibble())

# nrow(utr_source <- utr_source %>% distinct(target, gene_id)) # 31,654

# Join to targetscan out

df <- read_tsv(f)

df <- df %>% select(a_Gene_ID, miRNA_family_ID) %>%
  mutate(miRNA_family_ID = sapply(strsplit(miRNA_family_ID, "::"), `[`, 1) ) %>%
  separate(a_Gene_ID, into = c("target", "gene_id"), sep = ";") %>%
  dplyr::rename("query" = "miRNA_family_ID")

nrow(gene_features <- read_rds(paste0(wd, "gene_features.rds")))

df <- gene_features %>% right_join(df) 

paste_ids <- function(x) {
  x <- paste(x, sep = ';', collapse = ';')
}

df %>% 
  group_by(target) %>%
  summarise(across(query, .fns = paste_ids), n = n(), .groups = "drop_last") %>%
  view()



#. =====
str(three_prime_utr <- sapply(strsplit(x, " "), `[`, 1)) # three_prime_utr ids 


df %>% distinct(a_Gene_ID) # 861

df %>% distinct(miRNA_family_ID) # 262 mir:mir* 

# Q: how many arm-binding sites found?

df %>%
  mutate(MIRNA = "hairpin") %>%
  mutate(MIRNA = ifelse(grepl("mature", miRNA_family_ID), "mature", MIRNA)) %>%
  mutate(MIRNA = ifelse(grepl("star", miRNA_family_ID), "star", MIRNA)) %>%
  count(MIRNA)

# 

df <- df %>% 
  separate(a_Gene_ID, into = c("target","gene_id"), sep = ";") %>%
  separate(miRNA_family_ID, into = c("query", "coor"), sep = "::") 

# 
df <- df %>% 
  mutate(UTR = ifelse(target %in% three_prime_utr, "three_prime_utr", "five_prime_utr"))

df %>% count(UTR)

TARGETSCAN <- df %>%
  group_by(query) %>%
  count(sort = T)

df1 <- read_rds(paste0(wd, "/mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds"))

RNAHYBRID <- df1 %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  # filter(target %in% three_prime_utr) %>%
  group_by(query) %>%
  count(sort = T) 

RNAHYBRID %>% 
  left_join(TARGETSCAN, by = "query") %>% 
  rename("RNAHYBRID"="n.x", "TARGETSCAN"="n.y") %>%
  separate(query, into = c("Name", "arm"), sep = "[.]") %>%
  pivot_longer(cols = c("RNAHYBRID", "TARGETSCAN"), names_to = "Tool") %>%
  ggplot(aes(value, fill = Tool, color = Tool)) + 
  geom_histogram() +
  theme_bw(base_size = 12, base_family = "GillSans") +
  scale_x_continuous("Number of target", labels = scales::comma) +
  scale_y_continuous("Number of miRs", labels = scales::comma) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'right',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 

# SEARCH CONSISTENCY BETWEEN RNAhybrid AND TARGETSCAN

df1 <- df1 %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  filter(target %in% three_prime_utr)

str(q.RNAhybrid <- df1 %>% mutate(ID = paste(target, query, sep = "|")) %>% distinct(ID) %>% pull())
str(q.targetscan <- df %>% mutate(ID = paste(target, query, sep = "|")) %>% distinct(ID) %>% pull())


q.RNAhybrid[q.RNAhybrid %in% q.targetscan]
sum(q.RNAhybrid %in% q.targetscan)

nrow(df) # 13 418

nrow(df1) # 18 050

nrow(.DB <- df1 %>% inner_join(df, by = c("target", "query")))

view(.DB %>% drop_na(gene_id))

# OVERLAP w/ DB

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv"))

# 

DB %>% left_join(.DB)


# HEATMAP FROM ABUND.
df1 <- df1 %>% separate(query, into = c("Name", "arm"), sep = "[.]")

df1 <- df1 %>% drop_na(description) 

df1 <- df1 %>% distinct(gene_id, description, Name)

COUNTS %>% right_join(df1)

df1 %>% arrange(desc(description))
