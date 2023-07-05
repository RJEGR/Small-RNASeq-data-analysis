
# RICARDO GOMEZ-REYES  
# AFTER RUN 2_SRNA_LOCATION_DB_PREP.R

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

library(GenomicRanges)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

print(DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

view(DB)

# 
DB %>% count(biotype_best_rank, sort = T)

DB %>% count(SRNAtype, sort = T)

other_nc <- c("snRNA", "rRNA", "tRNA", "lncRNA", "srpRNA")

df <- DB %>% 
  mutate(width = nchar(MajorRNA)) %>%
  mutate(biotype_best_rank = ifelse(biotype_best_rank %in% other_nc, "Other ncRNA", biotype_best_rank)) %>%
  group_by(SRNAtype, biotype_best_rank, width) %>% 
  summarise(Reads = sum(Reads), n = n()) %>%
  group_by(SRNAtype) %>%
  mutate(reads_pct = Reads/sum(Reads), n_freq = n/sum(n))

df %>% summarise(sum(Reads), sum(n))

df %>% summarise(sum(reads_pct), sum(n_freq))


# (doi.org/10.1186/1471-2164-11-533): Roughly half of known miRNA genes are located within previously annotated protein-coding regions ("intragenic miRNAs"). A high-confidence set of predicted mRNA targets of intragenic miRNAs also shared many of these features with the host genes. Approximately 20% of intragenic miRNAs were predicted to target their host mRNA transcript.

df %>% filter(SRNAtype == "miR") %>% group_by(biotype_best_rank) %>% summarise(sum(n_freq))

df %>%
  ggplot(aes(x = SRNAtype, y = reads_pct, fill = biotype_best_rank)) +
  geom_col() +
  see::scale_fill_material() +
  scale_y_continuous("Freq. of reads",labels = scales::percent) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 
  # facet_grid( SRNAtype ~ .)

df %>%
  ggplot(aes(x = width, y = Reads, fill = biotype_best_rank)) +
  geom_col() +
  labs(x = "Read length (nt)") +
  see::scale_fill_material() +
  scale_y_continuous("Number of reads",labels = scales::comma) +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.minor = element_blank()) 
  # facet_grid(~ SRNAtype, scales = "free", space = "free")

# ANALYSIS OF CANONICAL OR ISOFORM (NEILSEN ET AL 2012) ?

DB %>% 
  drop_na(KnownRNAs) %>%
  ggplot(aes(MajorRNAReads/Reads)) +
  geom_histogram()

#

DB %>% 
  # drop_na(KnownRNAs) %>%
  ggplot(aes(FracTop, MajorRNAReads/Reads, color = Strand)) +
  facet_grid(~ SRNAtype) +
  geom_point()


# ANY CUTOFF?

DB %>% count(SRNAtype)

DB %>% 
  # filter()
  group_by(SRNAtype) %>%
  count(type, sort = T) %>%
  ggplot()
  
  
# 

.out %>% 
  # filter(Name %in% which_mirs) %>%
  count(Name, sort = T) 

# gggenes::example_dummies %>%
#   ggplot(aes(xmin = start, xmax = end, y = molecule,fill = gene)) +
#   geom_gene_arrow()

which_cluster <- "Cluster_2825"

cluster_bind <- assembly[assembly$ID == which_cluster] %>% as_tibble() %>% 
  rename("Name"="ID") %>%
  mutate(biotype = type, gene_id = DicerCall, transcript_id = MIRNA) %>% 
  select(any_of(names(.out)))


.out %>% 
  filter(Name == which_cluster) %>% 
  # filter(!type %in% c("mRNA", "gene", "ncRNA")) %>%
  rbind(cluster_bind) %>%
  ggplot(aes(xmin = start, xmax = end, y = 1, forward = T, fill = biotype)) +
  geom_gene_arrow() +
  facet_grid(type ~.)


  
  

