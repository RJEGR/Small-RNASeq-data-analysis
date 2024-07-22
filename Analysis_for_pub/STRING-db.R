# miRNA-mRNA cds vs STRING DB diamont-blast
# Processing results
# Considering using non-model orgs. result in low info per link
# So, let use a single mondel organism to retrieve protein-links

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

f <- list.files(dir, pattern = "protein.sequences.v12.0", full.names = T)


library(tidyverse)

read_outfmt6 <- function(f) {
  
  # seqid = transcript_id
  outfmt6.names <- c("transcript_id", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score")
  
  
  df <- read_tsv(f, col_names = F) 
  
  colnames(df) <- outfmt6.names
  
  df <- df %>% mutate(db = basename(f))
  
  return(df)
  
  
}

df <- do.call(rbind, lapply(f, read_outfmt6))

# df %>% ggplot(aes(identity, color = db)) + stat_ecdf()
df %>% ggplot(aes(x = identity, y = db)) + ggridges::geom_density_ridges()
# Join to species ----
# STRING_type: Core species are BLAST aligned all-against-all, periphery only against the core.
f <- list.files(dir, pattern = "STRING_DB_species.v12.0.txt", full.names = T)

taxon_df <- read_tsv(f) %>% mutate(taxon_id = as.character(`#taxon_id`))


df <- df %>% 
  mutate(taxon_id = sapply(strsplit(subject, "[.]"), `[`, 1)) %>% 
  left_join(taxon_df) 

df %>% count(official_name_NCBI, sort = T)

# df %>% filter(official_name_NCBI %in% "Lottia gigantea") %>% distinct(subject)

# Join to protein description ----

taxons <- df %>% distinct(taxon_id) %>% pull()

string_links <- function(x) {
  
  out_dir <- file.path(dir, "protein_links_full_v12_dir")
  
  system(paste0("mkdir -p ", out_dir))
  
  url <- "https://stringdb-downloads.org/download/protein.links.full.v12.0/"
  
  file <- paste0(x, ".protein.links.full.v12.0.txt.gz")
  
  file_url <- file.path(url, file)
  
  file_out <- file.path(out_dir, file)
  
  command <- paste0("wget -o ", file_out, " ", file_url)
  
  print(command)
  
  download.file(url = file_url, destfile = file_out)
  # file(command)
  

  }

# lapply(taxons, string_links)

links_dir <- file.path(dir, "protein_links_full_v12_dir")

# lapply(taxons[ !taxons %in% sapply(strsplit(list.files(links_dir), "[.]"), `[`, 1)], string_links)


string_info <- function(x) {
  
  
  # https://stringdb-downloads.org/download/protein.info.v12.0/225164.protein.info.v12.0.txt.gz
  
  out_dir <- file.path(dir, "protein_info_v12_dir")
  
  system(paste0("mkdir -p ", out_dir))
  
  url <- "https://stringdb-downloads.org/download/protein.info.v12.0/"
  
  file <- paste0(x, ".protein.info.v12.0.txt.gz")
  
  file_url <- file.path(url, file)
  
  file_out <- file.path(out_dir, file)
  
  # command <- paste0("wget -o ", file_out, " ", file_url)
  
  # print(command)
  
  download.file(url = file_url, destfile = file_out)
  # file(command)
  
  
}

info_dir <- file.path(dir, "protein_info_v12_dir")

lapply(taxons[ !taxons %in% sapply(strsplit(list.files(info_dir), "[.]"), `[`, 1)], string_info)

# lapply(taxons, string_info)



f <- list.files(file.path(dir, "protein_info_v12_dir"), pattern = "gz", full.names = T)


protein_info_df <- read_tsv(f, col_names = T) 

df <- df %>% left_join(protein_info_df, by = c("subject" = "#string_protein_id"))


# gene2tr -----

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

gene2tr <- read_rds(paste0(wd, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)


df <- df %>% left_join(gene2tr, by = "transcript_id")

names(df)

df <- df %>% distinct(gene_id, subject, preferred_name, taxon_id, official_name_NCBI, annotation)

write_tsv(df, file = file.path(dir, "gene2stringid_diamond_blastx.tsv"))

# Lot of data per file

# f <- list.files(file.path(dir, "protein_links_full_v12_dir"), pattern = "gz", full.names = T)

# links_df <- read_tsv(f[1], col_names = T) 
