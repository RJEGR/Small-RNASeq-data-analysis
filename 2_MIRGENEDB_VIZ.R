# RICARDO GOMEZ-REYES

# MICRORNA CONSERVATION
# Identificación de microRNAs conocidos en la etapa trocófora (pre-competente) y competente (veliger tardía) en el abulón rojo Haliotis rufescens a través de técnicas de secuenciación masiva.


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)


wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

dim(MIRGENEDB <- read_tsv(paste0(wd, "/MIRGENEDB_2.1.tsv")))

MIRGENEDB <- MIRGENEDB %>% select_at(vars(contains(c("MirGeneDB_ID", "Family","Node_of_origin_", "Seed"))))


# view(MIRGENEDB)


# 1) LOAD DATA ====
# COUNT THE NUMBER OF HOMOLOGY-BASED:SHORTSTACKS PREDICTED MIRS:

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out"

print(.DB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")))

# 1.1) SELECT ONLY KNOWN MIRS
.DB %>% filter(SRNAtype == "miR") %>% dplyr::count(SRNAtype) # 147

print(DB <- .DB %>% filter(SRNAtype == "miR") %>% drop_na(KnownRNAs))

DB %>% dplyr::count(SRNAtype) # MUST BE 69: 53TRUECONSERVED + 12HOMOLOGYBASED + 4 HOMOLOGYBASED PIR/MIR

# 1.2) SPLIT AND UNNEST MIR NAMES ====

DB <- DB %>%
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  dplyr::select(Name, KnownRNAs, MajorRNA) %>%
  unnest(KnownRNAs) %>%
  filter(!grepl("^piR-",KnownRNAs)) %>% # <- ALSO EXCLUDE OVERLAPED MIRS/PIRS
  separate(KnownRNAs, into = c("MirGeneDB_ID", "arm"), sep = "_") %>% 
  mutate(sp = sapply(strsplit(MirGeneDB_ID, "-"), `[`, 1)) %>%
  left_join(MIRGENEDB, by = "MirGeneDB_ID")

DB %>% distinct(Name) %>% nrow() # MUST BE 69 - 1: TRUECONSERVED BUT PIR HOMOLOGYBASED

# WRITE OUTPUT FOR FURTHER MERGING

write_tsv(DB, file = paste0(wd, "/SRNA2MIRGENEDB.tsv"))
  

# 1.3) SEED ANALYSIS: ====

DB %>% 
  mutate(SEED = substr(MajorRNA, 2,8) == Seed) %>%
  dplyr::count(SEED)

FALSE_SEED_DB <- DB %>% 
  mutate(SEED = substr(MajorRNA, 2,8) == Seed) %>%
  # mutate(SEED = grepl(paste(SEED, collapse = "|"), substr(MajorRNA, 2,8))) %>%
  filter(SEED != TRUE) %>% 
  left_join(.DB)

FALSE_SEED_DB %>% dplyr::count(Seed, sort = T)

# SOMETIMES, SEED WINDOW MATCH NOR BETWEEN 2-8 NUCL. POSITION
# EX:

WHICH_SEED <- c("UAAGGCA|UUGGUCC|CACAGCC|CCCUGUA|UUGUGAC|UUGCACU|UAGCACC|AUUGCUU")

FALSE_SEED_DB <- FALSE_SEED_DB %>% 
  filter(!grepl(WHICH_SEED, Seed))

FALSE_SEED_DB %>% dplyr::count(Name,Seed, Family, arm, sort = T)

# view(FALSE_SEED_DB)


# (https://doi.org/10.1093/nar/gkz885): Although highly conserved across vast distances of geologic time, seed sequences can and do change (11,23), expanding the functional repertoire of an ancestral seed sequence.

# THEN:

SEED <- MIRGENEDB %>% distinct(Seed) %>% pull()

MAJORRNA <- DB %>% distinct(MajorRNA) %>% pull()

sum(grepl(paste(SEED, collapse = "|"), MAJORRNA))

# CAUTION: FAMILIES ARE NOT NECESARY SINGLE SEED:

DB %>% distinct(Family, Seed) %>% filter(Family == "MIR-124")

DB %>% distinct(Family, Seed) %>% dplyr::count(Family, sort = T)


# 1.4) HOW MANY MIR FAMILIES FOUND:? ====

DB %>% distinct(Family) %>% nrow() # 49 PRINCIPAL FAMILIES

# ARE CLUSTERS-EXCLUSIVE FAMILIES FOUND 

DB %>% distinct(Name, Family) %>% dplyr::count(Family, sort = T)  

# SINGLE LOCUS FAMILIE:

str(DB %>% distinct(Name, Family) %>% dplyr::count(Family, sort = T) %>% filter(n == 1) %>% pull(Family))

# MULTI-LOCUS FAMILIES:

str(DB %>% distinct(Name, Family) %>% dplyr::count(Family, sort = T) %>% filter(n > 1) %>% pull(Family))


# 2) )DATA VIZ ====

# A) PLOT NUMBER OF FAMILIES PER CLUSTER:

LEFT_DB <- .DB %>% select(Name, type, biotype, WGCNA)

DB %>% distinct(Name, Family) %>% 
  # left_join(LEFT_DB) %>%
  # group_by(biotype) %>%
  dplyr::count(Family, sort = T) %>% 
  mutate(Family = fct_reorder(Family, n)) %>%
  ggplot(aes(y = Family, x = n)) +
  geom_col() +
  labs(x = "Number of clusters", y = "Family") +
  theme_classic(base_size = 10, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey', color = 'white')) -> ps

ggsave(ps, filename = 'FAMILIES_PER_CLUSTERS.png', 
  path = wd, width = 4, height = 8, device = png, dpi = 300)


# B) PLOT NUMBER SPECIES PER FAMILIE
# HOW CONSERVED (BY FAMILY)?

CONSERVED_TO_HALIOTIS <- DB %>% 
  mutate(Node = `Node_of_origin_(family)`) %>%
  distinct(Name, Node) %>%
  dplyr::count(Node, sort = T) %>%
  rename("Haliotis"="n")

MIRGENEDB %>% 
  mutate(Node = `Node_of_origin_(family)`) %>%
  dplyr::count(Node, sort = T) %>%
  rename("DB"="n") %>% 
  right_join(CONSERVED_TO_HALIOTIS) %>% 
  view()


# HOW MANY CLADES?

DB %>% 
  distinct(Family, `Node_of_origin_(family)`) %>%
  mutate(Node = `Node_of_origin_(family)`) %>%
  dplyr::count(Node, sort = T)


# DB %>% filter(sp == "Lgi") %>% distinct(Family) 
# DB %>% filter(`Node_of_origin_(family)` == "Lophotrochozoa") %>% distinct(Family) 


NodeLev <- DB %>% 
  mutate(Node = `Node_of_origin_(family)`) %>%
  group_by(Node) %>%
  distinct(sp) %>%
  dplyr::count(sort = T) %>% pull(Node)


# HOW MANY SPECIES?
DB %>% distinct(sp) %>% nrow() # 71

LEFT_DB <- DB %>% 
  distinct(sp, `Node_of_origin_(family)`) %>% 
  group_by(sp) %>%
  # count(`Node_of_origin_(family)`, sort = T)
  sample_n(1)

# INSTEAD OF CLUSTERS TRY FAMILIES

DB %>% 
  distinct(Family, sp) %>% 
  count(sp) %>%
  left_join(LEFT_DB) %>%
  mutate(sp = fct_reorder(sp, n)) %>%
  # mutate(`Node_of_origin_(family)` = factor(`Node_of_origin_(family)`, levels = NodeLev)) %>%
  ggplot(aes(y = sp, x = n)) +
  facet_grid(`Node_of_origin_(family)` ~., scales = "free_y", space = "free_y") +
  geom_col() +
  labs(x = "Number of families", y = "Organism") +
  theme_classic(base_size = 10, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey', color = 'white')) -> ps



ggsave(ps, filename = 'SPECIES_PER_FAMILIES.png', 
  path = wd, width = 4, height = 8, device = png, dpi = 300)
