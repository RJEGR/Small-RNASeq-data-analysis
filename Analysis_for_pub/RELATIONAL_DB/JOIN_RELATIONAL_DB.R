
# 0) LOAD gene2tr DATABASE
# 1) LOAD AND PARSE SMP RESULTS FROM BLAST miRNA:mRNA transcripts to SMP DB
# 2) LOAD EGGMAPPER RESULTS FROM miRNA:mRNA transcripts
# 3) LOAD STRING AND PARSE
# 4) LOAD miRNA:mRNA  TARGET DB 
# 4) JOIN DB
# LOAD DEGS COLS

library(tidyverse)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

out_dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

# 0
dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

genedescr <- read_rds(paste0(dir, "genome_features.rds"))[[2]]

gene2tr <- read_rds(paste0(dir, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)

# 1) Shell Matrix Proteins ----

dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = ".outfmt6", full.names = T)

# selecting only biomineralization
f <- f[grepl("Combined_Gastropoda|shell_matrix", f)]

library(tidyverse)

paste_col <- function(x) { 
    x <- x[!is.na(x)] 
    x <- unique(sort(x))
    x <- paste(x, sep = ';', collapse = ';') 
  
  return(x)

}

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

SMP_DB <- do.call(rbind, lapply(f, read_outfmt6))

SMP_DB <- SMP_DB %>% distinct(transcript_id, subject) %>% left_join(gene2tr)

SMP_DB <- SMP_DB %>% distinct(gene_id, subject) %>% dplyr::rename("SMPID" = "subject")

SMP_DB <- SMP_DB %>%
  group_by(gene_id) %>%
  summarise(across(SMPID, .fns = paste_col))

# 2) MAPPER ----
dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = "eggnog_mapper.emapper.annotations", full.names = T)

MAPPER_DB <- read_tsv(f, comment = "##")

eggNOG_cols <- c("gene_id", "Preferred_name","Description", "GO","PFAMs","BRITE", "CAZy", "COG_category")

MAPPER_DB <- gene2tr %>%
  right_join(MAPPER_DB, by = c("transcript_id" = "#query" )) %>%
  select(- transcript_id, -evalue, -score) %>%
  select_at(vars(contains(eggNOG_cols), starts_with("KEGG"))) %>%
  distinct()

MAPPER_DB %>% distinct(gene_id)
# Bind to the COG_category

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/"

NOG.col <- read_rds(paste0(dir, '/cogs.rds')) %>%
  dplyr::rename("COG_name" = "name")

MAPPER_DB <- MAPPER_DB %>% 
  left_join(NOG.col, by = c("COG_category" = "code"))
 
# 3) STRING----

dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = "protein.sequences.v12.0", full.names = T)

f <- f[grepl("9606", basename(f))] # Only human ID

STRING_DB <- do.call(rbind, lapply(f, read_outfmt6))

f <- list.files(dir, pattern = "STRING_DB_species.v12.0.txt", full.names = T)

taxon_df <- read_tsv(f) %>% mutate(taxon_id = as.character(`#taxon_id`))

STRING_DB <- STRING_DB %>% 
  mutate(taxon_id = sapply(strsplit(subject, "[.]"), `[`, 1)) %>% 
  left_join(taxon_df) 

f <- list.files(file.path(dir, "protein_info_v12_dir"), pattern = "gz", full.names = T)

f <- f[grepl("9606", basename(f))]

protein_info_df <- read_tsv(f, col_names = T) 

STRING_DB <- STRING_DB %>% left_join(protein_info_df, 
  by = c("subject" = "#string_protein_id"))

STRING_DB <- STRING_DB %>% left_join(gene2tr, by = "transcript_id")

STRING_DB <- STRING_DB %>% distinct(gene_id, preferred_name) %>% 
  dplyr::rename("STRINGID" = "preferred_name")

STRING_DB <- STRING_DB %>%
  group_by(gene_id) %>%
  summarise(across(STRINGID, .fns = paste_col))

# 4) LOAD miRNA:mRNA and MIR-LOCI ----
dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(.TARGETDB <- read_rds(paste0(dir,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))

is_na <- function(x) ifelse(is.na(x), 0, x)

keep_expressed <- .TARGETDB %>% 
  dplyr::select(starts_with("SRR")) %>% 
  mutate(across(where(is.double), ~is_na(.))) %>% 
  rowSums()

keep_expressed <- keep_expressed > 1

.TARGETDB <- .TARGETDB[keep_expressed,]

TARGETDB <- .TARGETDB %>% select(query, gene_id, description)

TARGETDB %>% distinct(gene_id, description) %>%
  rename("Description_cds" = "description") %>%
  left_join(MAPPER_DB) %>%
  rename("Description_mappeper" = "Description") %>%
  distinct(gene_id, Description_cds, Description_mappeper, GOs) %>%
  write_tsv(file.path(out_dir, "Gene_ontologies_DB.tsv"))

  
# MIR LOCI

dir <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.LOCIDB <- read_rds(paste0(dir, "RNA_LOCATION_MIR_DB.rds")))

MajorRNA2Name <- .LOCIDB %>% distinct(MajorRNA, Name) %>%   dplyr::rename("query" = "Name")

# Replace ClusterID w/ MajorRNA

TARGETDB <- TARGETDB %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  # distinct(gene_id, query) %>%
  mutate(query = gsub(".mature", "", query)) %>%
  left_join(MajorRNA2Name) %>% select(-query) %>%
  select(gene_id, MajorRNA, description)

# 5) JOIN DB -----
# Example:
TARGETDB %>% left_join(SMP_DB) %>% drop_na(SMPID)

DB <- TARGETDB %>% left_join(SMP_DB) %>% left_join(STRING_DB) 

nrow(DB %>% distinct(gene_id)) # MUST MATCH 142 ASTRINGENT GENES

eggNOG_cols <- c("COG_category", "COG_name")

DB <- MAPPER_DB %>% 
  filter(COG_category != "S") %>% # Function unknown
  select_at(vars(contains(c("gene_id", eggNOG_cols)))) %>%
  distinct() %>% drop_na(COG_name) %>% 
  right_join(DB)

nrow(DB %>% distinct(gene_id))

DB %>% distinct(gene_id, MajorRNA)


# LOAD MICRORNA DEGS -----

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

DEGSDB <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1)

RNA2ID <- DEGSDB %>% distinct(MajorRNA, MirGeneDB_ID) %>%
  dplyr::rename("MajorRNAID" = "MirGeneDB_ID")

DB <- DEGSDB %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast)) %>%
  group_by(MajorRNA) %>%
  summarise(across(Contrast, .fns = paste_col)) %>%
  left_join(RNA2ID) %>%
  right_join(DB) 
  
  
write_tsv(DB, file = file.path(out_dir, "LONGER_RELATIONAL_DB.tsv"))

# esquisser::esquisser(DB)



# EXIT ----

# 6) QUERY DEGS ----

DB %>% filter(STRINGID %in% which_targets)

graph

# 7) UPGRADE GOs. from genes ----
# Replace GOs from Swissprot by eggnogg
# See n GOs are higher in eggnog res.

TARGETDB %>% distinct(gene_id, MajorRNA)

DB %>% distinct(gene_id, MajorRNA) %>%
  left_join(distinct(MAPPER_DB, gene_id, GOs)) %>%
  distinct(MajorRNA, GOs) %>%
  

# 7.1) If want to bind: Gene Ontology processing ----

.TARGETDB %>% 
  drop_na(GO.ID) %>%
  distinct(gene_id, GO.ID) %>%
  mutate(GOs = strsplit(GO.ID, ";")) %>%
  unnest(GOs)

GODB <- MAPPER_DB %>% distinct(gene_id, GOs) %>%
  filter(GOs != "-") %>% 
  mutate(GOs = strsplit(GOs, ",")) %>%
  unnest(GOs) %>%
  group_by(gene_id) %>%
  summarise(across(GOs, .fns = paste_col))

GODB <- GOSWISS %>% select(-n_swiss) %>%
  rbind(GODB) %>%
  mutate(GOs = strsplit(GOs, ";")) %>%
  unnest(GOs) %>%
  group_by(gene_id) %>%
  summarise(across(GOs, .fns = paste_col))
  
DB %>% left_join(GODB) %>% view()
