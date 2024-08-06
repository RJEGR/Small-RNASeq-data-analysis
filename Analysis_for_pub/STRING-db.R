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

lapply(taxons[ !taxons %in% sapply(strsplit(list.files(links_dir), "[.]"), `[`, 1)], string_links)


string_info <- function(x) {
  
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

# Lot of data per file when using best hit (all sequences source)
# using only human (Model, 9606) and c. elegans (Model of larval development, 6239).
dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

graph <- read_rds(file.path(dir, "protein_links_full_v12.rds"))

f <- list.files(file.path(dir, "protein_links_full_v12_dir"), pattern = "gz", full.names = T)

f <- f[grepl("9606", basename(f))] # |6239

links_df <- read_delim(f, col_names = T, delim = " ")

query_nodes <- df %>% distinct(subject) %>% pull() %>% sort()

str(protein1 <- links_df %>% pull(protein1) %>% sort())
str(protein2 <- links_df %>% pull(protein2) %>% sort())

sum(protein1 <- protein1 %in% query_nodes)
sum(protein2 <- protein2 %in% query_nodes)

identical(protein1, protein2)

links_df_ <- links_df[protein1,]

# to hard to run
# nodes_list <- paste(query_nodes, collapse = "|")
# links_ <- c("protein1", "protein2")
# links_df_ <- links_df %>% filter_at(vars(all_of(links_)), any_vars(grepl(nodes_list, .)))

Node1 <- links_df_ %>% distinct(protein1) %>% 
  dplyr::rename("Node" = "protein1")

Node2 <- links_df_ %>% distinct(protein2) %>% 
  dplyr::rename("Node" = "protein2")

nrow(Nodes <- rbind(Node1, Node2))

Nodes <- Nodes %>% left_join(df, by = c("Node" = "subject"))

# subset for viz
which_targets <- c("RNF166", "CENPM", "CDK10","H2AZ2", "C1GALT1", "TFB1M","ABCG2","MAP11")

visNetwork::visNetwork(nodes = Nodes, edges = links_df_)

library(igraph)
library(tidygraph)
library(ggraph)

graph <- tidygraph::tbl_graph(nodes = Nodes, edges = links_df_, directed = T)

# graph <- graph %>% activate("nodes") %>% filter(Node %in% query_nodes)

write_rds(list(graph, df, links_df_), file = file.path(dir, "protein_links_full_v12.rds"))


graph <- graph %>% activate("nodes") %>% 
  mutate(Col = ifelse(preferred_name %in% which_targets, "Target under OA", ""),
         label = ifelse(preferred_name %in% which_targets, preferred_name, ""))
  
centrality_df <- graph %>% activate("nodes") %>% 
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
    # membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
    pageRank = page_rank(.)$vector) %>% as_tibble() %>% 
  select(Node, preferred_name, degree, betweenness, pageRank) %>% 
  arrange(degree) 

centrality_df <- centrality_df %>% filter(degree > 0) 

layout <- graph %>% activate("nodes") %>% 
  # filter(taxon_id == "9606") %>%
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
        # membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
        pageRank = page_rank(.)$vector) %>%
  filter(degree > 0) %>%
  create_layout( layout = 'kk')
  # create_layout( layout = 'linear', circular = TRUE)

# if linear
# ggraph(layout) +
#   geom_edge_arc(aes(alpha = after_stat(index)), strength = 0.2) + 
#   geom_node_text(aes(label = preferred_name, angle = node_angle(x, y), hjust = x < 0, color = Col)) + 
#   coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))

ggraph(layout) +
  geom_edge_link(aes(alpha = combined_score)) +
  geom_node_point(aes(color = Col, size = degree)) +
  geom_node_text(aes(label = preferred_name), repel = TRUE, family = "GillSans")
  # facet_nodes(~ official_name_NCBI, scales = "free")


# try
gr <- create_notable('Meredith') %>%
  mutate(class = sample(c('Class 1', 'Class 2', 'Class 3'), n(), replace = TRUE))

ggraph(gr, 'linear', circular = TRUE) +
  geom_edge_arc(aes(alpha = after_stat(index)), strength = 0.2) + 
  geom_node_text(aes(label = class, angle = node_angle(x, y), hjust = x < 0)) + 
  coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2))
