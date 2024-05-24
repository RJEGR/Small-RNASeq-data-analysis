
# LOAD REGULATORY FUNCTION DB
# LOAD Results.txt (shortstacks) 
# bind MAjorRNA column to DB using Cluster_ Name from Results.txt
# LOAD Upgraded Names from SEQUENCES_MERGED
# Bind target db (description, gene_id, etc.) to SEQUENCES_MERGED

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(TARGETDB <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv")))

TARGETDB <- TARGETDB %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  separate(query, into = c("query", "arm"), sep = "[.]") %>% 
  filter(arm == "mature") %>%
  filter(predicted == "BOTH") %>%
  select(query, gene_id, GO.ID, description, type, biotype, target) %>% rename("Name" = "query")

TARGETDB %>% count(Name)
TARGETDB %>% count(biotype)
TARGETDB %>% count(gene_id)

TARGETDB %>% drop_na(GO.ID) %>% count(gene_id)
TARGETDB %>% count(gene_id)

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

TARGETDB <- read_tsv(paste0(wd, "Results.txt")) %>%
  select(Name, MajorRNA) %>%
  right_join(TARGETDB, by = "Name") %>% select(-Name)


wd <-  "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

RES.P <- read_tsv(list.files(path = wd, pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>% 
  filter( padj < 0.05  & abs(log2FoldChange) > 1)

DESEQRES <- RES.P %>% distinct(MajorRNA, Name) # %>% rename("NewName" = "Name")

TARGETDB <- DESEQRES %>% right_join(TARGETDB, by = "MajorRNA")

# dds <- read_rds(paste0(wd, "/SEQUENCES_MERGED_DDS_DESEQ2.rds"))

TARGETDB %>% count(Name)

write_tsv(TARGETDB, paste0(path_out, "SEQUENCES_MERGED_MIRNA_TARGET_DB.tsv"))

# 2)

WHICH_CONT <- c("CONTRAST_A", "CONTRAST_B")

# Use description to 

RES.P <- TARGETDB %>% distinct(MajorRNA, gene_id, description) %>%
  # rename("target_id" = "gene_id")
  group_by(MajorRNA) %>%
  summarise(across(gene_id, .fns = paste_col), n = n()) %>%
  right_join(RES.P, by = "MajorRNA")

DESCRIPTION <- RES.P %>% 
  # filter by something
  # filter(CONTRAST %in% WHICH_CONT) %>%
  select(gene_id, Name) %>%
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id) %>%
  distinct(gene_id) %>%
  left_join(TARGETDB %>% distinct(gene_id, description), by = "gene_id")


# MIRNA LOCI REGULATING ITS HOST? (ZERO)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

LOCUSDB <- read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>% 
  dplyr::filter(SRNAtype == "miR") %>%
  mutate(gene_id = ifelse(is.na(gene_id), Locus, gene_id)) %>%
  dplyr::select(MajorRNA, gene_id, biotype_best_rank, type, biotype, transcript_id)
  # distinct(MajorRNA)

LOCUSDB %>% count(biotype_best_rank)

gene_id_target <- RES.P %>% 
  select(gene_id) %>%
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id) %>% arrange(gene_id) %>% pull(gene_id)

gene_id_locus <- LOCUSDB %>% 
  filter(biotype_best_rank == "Intragenic") %>%
  arrange(gene_id) %>% pull(gene_id)

sum(gene_id_target %in% gene_id_locus)

# IS THE TARGER USUALLY EXPRESED DURING LARVAL DEV.

wd <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024"
dim(MIR_COUNT <- read_rds(file.path(wd, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds"))) #  must match 117 rows x 12 cols

MIR_COUNT <- rowSums(MIR_COUNT) %>% as_tibble(rownames = 'MajorRNA') %>% rename("mir_exp" = "value")

# view(MIR_COUNT)

MIRNA2TARGET <- RES.P %>% 
  # filter by something
  # filter(CONTRAST %in% WHICH_CONT) %>%
  select(gene_id, MajorRNA) %>%
  mutate(gene_id = strsplit(gene_id, ";")) %>%
  unnest(gene_id) %>%
  left_join(MIR_COUNT, by = "MajorRNA")

# sanity check

MIRNA2TARGET %>% count(gene_id)

n_mirs <- function(x) {length(unique(x))}

MIRNA2TARGET <- MIRNA2TARGET %>%
  group_by(gene_id) %>%
  summarise(n_mirs = n_mirs(MajorRNA), across(MajorRNA, .fns = paste_col), mir_exp = sum(mir_exp)) 

# LEER CONTEOS DE GENE-EXPRESSION RNA-SEQ Y SUMAR INDEPENDIENTE DE LA ETAPA DE DESARROLLO

GTF <- "transcripts.gtf"

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"

f <- list.files(path = wd, pattern = "METADATA_RNASEQ", full.names = T)

.colData <- read_csv(f)

f <- list.files(path = wd, pattern = GTF, full.names = T)

GTF <- rtracklayer::import(f)

print(GTF2DF <- GTF %>% as_tibble() %>% 
    distinct(gene_id, transcript_id, ref_gene_id) %>% 
    filter(!is.na(ref_gene_id))) # A tibble: 70,001 × 3

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/QUANTIFICATION/"

f <- list.files(path = wd, pattern = "gene_count_matrix.csv", full.names = T)

dim(.COUNT <- read_csv(f))

keep_sam <- .colData %>% filter(grepl("dpf", CONTRAST_B)) %>% pull(LIBRARY_ID)

dim(.COUNT <- .COUNT %>% select(all_of(c("gene_id", keep_sam))))

.COUNT %>% distinct(gene_id) %>% nrow() # 58592 genes assembled

any(GTF2DF$gene_id %in% .COUNT$gene_id)

GTF2DF <- GTF2DF[GTF2DF$gene_id %in% .COUNT$gene_id,]

row_names <- .COUNT$gene_id

head(.COUNT <- data.frame(.COUNT, row.names = row_names) %>% select(-gene_id))

TARGET_COUNT <- rowSums(.COUNT) %>% as_tibble(rownames = 'gene_id') %>% rename("target_exp" = "value")

TARGET_COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>%
  right_join(TARGET_COUNT) %>%
  rename("assembly_id" = "gene_id", "gene_id" = "ref_gene_id") %>%
  # select(-assembled_id) %>%
  distinct()

TARGET_COUNT %>% count(gene_id, sort = T)

MIRNA2TARGET <- TARGET_COUNT %>% right_join(MIRNA2TARGET, by = "gene_id") %>% right_join(DESCRIPTION)
 
write_tsv(MIRNA2TARGET, paste0(path_out, "SEQUENCES_MERGED_MIRNA2TARGET_EXP_DB.tsv"))

query.targets <- MIRNA2TARGET %>% distinct(gene_id) %>% pull()

COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>%
  right_join(read_csv(f)) %>%
  rename("assembly_id" = "gene_id", "gene_id" = "ref_gene_id") %>%
  filter(gene_id %in% query.targets)

nrow(COUNT) # 155 from the 168 targets usually are expressed during larval dev.

write_tsv(COUNT, paste0(path_out, "REFBASED_MODE_COUNT.tsv"))

# heatmap(COUNT %>% 
#   dplyr::select(contains(c("SRR"))) %>%
#   as("matrix"))

cor.mat <- COUNT %>% 
  # mutate_at(vars(contains(c("SRR"))), z_scores) %>%
  left_join(MIRNA2TARGET) %>%
  drop_na(target_exp) %>%
  filter(target_exp > 0) %>%
  # mutate(target_exp = log10(target_exp), mir_exp = log10(mir_exp/n_mirs)) %>%
  mutate(target_exp = log2(target_exp), mir_exp = log2(mir_exp/n_mirs)) %>%
  mutate(ratio = mir_exp-target_exp) %>%
  dplyr::select(contains(c("n_mirs", "mir_exp", "target_exp", "ratio"))) %>%  # "SRR", 
  rstatix::cor_mat(method = "spearman")

cor.mat %>% cor_get_pval() %>% as_tibble(rownames = "rowname") %>%
  pivot_longer(-rowname, values_to = "cor_pval") %>% arrange(rowname, name)  %>% View()

cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)

MIRNA2TARGET %>%
  drop_na(target_exp) %>%
  filter(target_exp > 0) %>%
  arrange(desc(n_mirs)) %>%
  # mutate(target_exp = log10(target_exp), mir_exp = log10(mir_exp)) %>% # or mir_exp = log10(mir_exp/n_mirs)
  mutate(target_exp = log2(target_exp), mir_exp = log2(mir_exp/n_mirs)) %>%
  mutate(ratio = mir_exp-target_exp) %>%
  # mutate(ratio = mir_exp/target_exp) %>%
  ggplot(aes(target_exp, ratio)) + geom_point(alpha = 0.5, aes(size = n_mirs), shape = 21) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  labs(x = "Target expression (log2)", y = "microRNA fold-change")


# mutual inf analysis and anticorrelation ()

dataset <- COUNT %>% 
  # mutate_at(vars(contains(c("SRR"))), z_scores) %>%
  left_join(MIRNA2TARGET) %>%
  drop_na(target_exp) %>%
  filter(target_exp > 0) %>%
  mutate(target_exp = log2(target_exp), mir_exp = log2(mir_exp/n_mirs)) %>%
  mutate(ratio = mir_exp-target_exp) %>%
  dplyr::select(contains(c("n_mirs", "mir_exp", "target_exp", "ratio"))) # %>%  # "SRR", 

names(dataset)[names(dataset) %in% .colData$LIBRARY_ID] <- .colData$CONTRAST_B

corm <- cor(dataset, method = "spearman") # %>% cor_reorder()

heatmap(corm)

corm <- corm %>% as_tibble(rownames = "rowname") %>%
  pivot_longer(-rowname, values_to = "cor") %>% arrange(rowname, name)

library(minet)

# Build Mutual Information Matrix

mim <- minet::build.mim(dataset, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset)))

heatmap(mim)

# Context Likelihood or Relatedness Network
clr(mim)
net<-mrnet(mim)

heatmap(net)

res<-norm(net)

table <- validate(res, net)

df <- mim %>% as_tibble(rownames = "rowname") %>% 
  pivot_longer(-rowname, values_to = "mim") %>% 
  arrange(rowname, name) %>% left_join(corm)
  
View(df)  

df %>% 
  filter(rowname == "ratio") %>%
  ggplot(aes(mim, cor)) + 
  geom_point(size = 0.5) + geom_text(aes(label = name)) + 
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  # ggpubr::stat_cor(method = "spearman", cor.coef.name = "R", p.accuracy = 0.001) +
  facet_wrap(~ rowname)

view(df)

