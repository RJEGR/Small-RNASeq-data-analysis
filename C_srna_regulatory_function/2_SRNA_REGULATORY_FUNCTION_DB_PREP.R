# RICARDO GOMEZ REYES
# CREAR UN CODIGO LIMPIO EN EL QUE OCURRA LO SIGUIENTE:
# 1) LEA LOS ARCHIVOS DE RNAHYBRID Y TARGETSCAN
# 2) CURE LOS RESULTADOS A UN FORMATO COMPATIBLE PARA HACER BIND DE AMBOS OBJETOS
# 3) CREE UNA COLUMNA NUEVA INDIQUE SI LA PREDICCION FUE RESUELTA EN AMBAS HERRAMIENTAS O ALGUNA
# 4) ALMACENE UN ARCHIVO DE NOMBRE SRNA_FUNCTION_PREDICTED.TSV PARA SUBSECUENTES ANALISIS

# LAS COLUMNAS CONTENIDAS DEBEN SER: 
# seqnames, 
# gene_coords (start:end:strand), 
# gene_id (related to UTR) <- key value for step 3_ 
# target_id (Ids from UTR.fasta used during the target analysis)
# query_id (Names from shortstacks DB)
# n (# srnas, if unique target n == 1, else UTR is binding site by multiple mirs)
# predicted_at (targetscan, rnahybrid, or both)

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

# 1) ====

pattern <- "mir_vs_utr_rmdup_RNAhybrid.out.psig.tsv"

f <- list.files(path = wd, pattern = pattern, full.names = T)

# AFTER CLEAN RNAHYBRID pval < 0.05, JUST LOAD .out.psig.tsv file:

RNAHYBRID <- read_tsv(f)

# TARGETSCAN INCLUDE ONLY 3 UTR FROM psig RNAHybrid detection

f <- "mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out"

f <- list.files(path = wd, pattern = f, full.names = T)

TARGETSCAN <- read_tsv(f)

# BIND TO TARGETSCAN, UPGRADE VERSION OF RESIDUAL 17 MIRS

f <- "mir_ts_vs_three_prime_utr_rmdup_ts_targetscan.out"

f <- list.files(path = wd, pattern = f, full.names = T)

TARGETSCAN <- read_tsv(f) %>%
  mutate(miRNA_family_ID = paste0(miRNA_family_ID, ".mature")) %>%
  rbind(TARGETSCAN, .)


nrow(TARGETSCAN %>% distinct(miRNA_family_ID) %>% filter(grepl("mature", miRNA_family_ID))) # MUST BE 147


# 2) ====
# 2.1) PREPARE UTR INFO ====

# from gene feature to UTR (flat information)

utr_f <- list.files(path = wd, full.names = T, pattern = "three_prime_utr.ids")

str(x <- read_lines(utr_f)) # 54 432 <-- i.e the N sequences from three_prime_utr.fa

str(gene_id <- sapply(strsplit(x, " "), `[`, 2)) # LOC*

str(target <- sapply(strsplit(x, " "), `[`, 1)) # three_prime_utr ids 

nrow(utr_source <- data.frame(target, gene_id) %>% as_tibble()) # 54 432

nrow(utr_source <- utr_source %>% distinct(target, gene_id)) # 31,654


# 2.2) =====
utr <- sapply(strsplit(x, " "), `[`, 1)

RNAHYBRID <- RNAHYBRID %>%
  mutate(query = sapply(strsplit(query, "::"), `[`, 1) ) %>%
  select(target, query) %>%
  filter(target %in% utr) %>% # <- keep only 3'utr
  mutate(predicted = "RNAHYBRID") 

# 2.3)


TARGETSCAN <- TARGETSCAN %>% select(a_Gene_ID, miRNA_family_ID) %>%
  mutate(miRNA_family_ID = sapply(strsplit(miRNA_family_ID, "::"), `[`, 1) ) %>%
  separate(a_Gene_ID, into = c("target", "gene_id"), sep = ";") %>%
  dplyr::rename("query" = "miRNA_family_ID") %>%
  select(target, query) %>%
  mutate(predicted = "TARGETSCAN")

identical(names(TARGETSCAN), names(TARGETSCAN))

# 3) ====

paste_ids <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';') 
  }

which_tools <- function(x) { 
  x <- x[!is.na(x)] 
  n <- length(unique(x))
  x <- unique(x)
  
  if(n > 1) {
    x <- "BOTH"
  } else
  x <- paste(x, sep = '|', collapse = '|') }


out <- rbind(RNAHYBRID, TARGETSCAN) %>%
  # sample_n(100) %>%
  group_by(target) %>%
  summarise(
    across(query, .fns = paste_ids), 
    across(predicted, .fns = which_tools), 
    n_rnas = n(),
    .groups = "drop_last")

out %>% count(predicted)

nrow(out) # 9565 different 3' utrs predicted to be target by srna

# ES NECESESARIO CARGARA LAS ANOTACIONES  A NIVEL GENOMA Y TRANSCRIPTOMA:

nrow(out <- out %>% left_join(utr_source) %>% arrange(desc(n_rnas)))

# out <- out %>% select(target, query, n_rnas, predicted)

write_rds(out, file = paste0(wd, "SRNA_FUNCTION_PREDICTED.rds"))


