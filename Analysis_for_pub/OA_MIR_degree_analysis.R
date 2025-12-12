
library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

RES <- read_tsv(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>%
  mutate(Contrast = ifelse(sign(log2FoldChange) == -1, sampleB, sampleA)) %>%
  mutate(Contrast = paste0(CONTRAST,"_",Contrast))

# RES.P <- RES %>% filter( padj < 0.05  & abs(log2FoldChange) > 1)



QUERIES <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B")) %>%
  # mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  distinct(MajorRNA, MirGeneDB_ID) %>%
  pull(MajorRNA, name = MirGeneDB_ID)

# QUERIES

DF1 <- DB %>% 
  # mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>% 
  distinct(STRINGID, gene_id, MajorRNA, MajorRNAID) %>%
  count(STRINGID, gene_id, sort = T) %>%
  dplyr::rename("MIR_degree" = "n")
  


paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

DF2 <- DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  distinct(gene_id, MajorRNAID) %>%
  mutate(Contrast = ifelse(MajorRNAID %in% names(QUERIES[1:3]), "24 hpf", "110 hpf")) %>%
  group_by(gene_id, Contrast) %>%
  summarise(across(MajorRNAID, .fns = paste_col), OA_MIR_degree = n()) %>%
  arrange(desc(OA_MIR_degree)) %>%
  mutate(OA_MIR_degree)
  # dplyr::rename("preferred_name" = "STRINGID")


NodeDF <- DF2 %>% 
  right_join(DF1 , by = "gene_id") %>%   
  mutate(OA_MIR_degree = ifelse(is.na(OA_MIR_degree), 0, OA_MIR_degree))

NodeDF <- DB %>% 
  # mutate(COG_name = ifelse(is.na(COG_name), "Unknown", STRINGID)) %>%
  mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
  distinct(gene_id, COG_name) %>%
  mutate(COG_name = factor(COG_name, levels = sort(decreasing = T, unique(COG_name)))) %>%
  right_join(NodeDF, by = "gene_id")

# Plot something intersting
# by protein

NodeDF %>%
  mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>% 
  mutate(Frac = MIR_degree) %>%
  filter(Frac > 0) %>% 
  group_by(Contrast) %>%
  arrange(desc(Frac), .by_group = T) %>% 
  mutate(Label = STRINGID, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
        levels = rev(paste(Label, row_number, sep = "__")))) %>%
  # mutate(Label = factor(Label, levels = unique(Label))) %>%
  ggplot(aes(y = Label, x = Frac, color = COG_name)) +
  # facet_wrap(~ Contrast, nrow = 3, ncol = 1, scales = "free_y") + # nrow = 2, ncol = 1, 
  facet_grid( Contrast ~ ., space = "free_y", scales = "free_y") +
  geom_point() +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) 

# by functional 

NodeDF %>%
  group_by(Contrast, COG_name) %>%
  summarise(MIR_degree = sum(MIR_degree), OA_MIR_degree = sum(OA_MIR_degree)) %>%
  mutate(Frac = MIR_degree) %>%
  filter(Frac > 0) %>% 
  group_by(Contrast) %>%
  arrange(desc(Frac), .by_group = T) %>% 
  mutate(Label = COG_name, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  # mutate(Label = factor(Label, levels = unique(Label))) %>%
  ggplot(aes(y = Label, x = Frac, color = Contrast, group = Contrast)) +
  # facet_wrap(~ Contrast, nrow = 3, ncol = 1, scales = "free_y") + # nrow = 2, ncol = 1, 
  # facet_grid( Contrast ~ ., space = "free_y", scales = "free_y") +
  geom_point() +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x))
  
# Visualice network

dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

print(stringg <- read_rds(file.path(dir, "protein_links_full_v12.rds")))

# Processing edge evidence
library(tidygraph)
library(igraph)

edges <- stringg %>% activate("edges") %>% as_tibble()

evidence <- data.frame(edges$fusion > 0,
  edges$neighborhood > 0,
  edges$cooccurence > 0,
  edges$experiments > 0,
  edges$textmining > 0,
  edges$database > 0,
  edges$coexpression > 0)

stringg <- stringg %>% activate("edges") %>% 
  mutate(Evidence_size = rowSums(evidence)/ncol(evidence))

stringg <- stringg %>% activate("nodes") %>% 
  filter(centrality_degree() > 0)

stringg <- stringg %>% activate("nodes") %>% 
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
    pageRank = page_rank(.)$vector) 

stringg <- stringg %>% activate("nodes") %>% left_join(NodeDF, by = "gene_id") 

# stringg <- stringg %>% activate("nodes") %>% as_tibble()

# Get edges connected to the specified node ----

QUERY <- NodeDF %>% drop_na(MajorRNAID) %>% view()

str(nodes_id <- which(V(stringg)$preferred_name %in% QUERY))

common_connected_edges <- incident(stringg, nodes_id)

# Get all neighbors of the specified node
common_vertices <- neighbors(stringg, nodes_id)
print("Connected Vertices:")
print(common_vertices)

# walk_paths <- all_shortest_paths(stringg, from = common_vertices)

# data.frame(path = unique(unlist(walk_paths$res)))

common_vertices <- V(stringg)$preferred_name[c(common_vertices)]

all_neighbors_shortest_paths <- function(graph, v) {
  
  # connected_edges <- incident(stringg, node_id)
  Vname <- v
  
  v <- which(V(graph)$preferred_name %in% v)
  
  # Get all neighbors of the specified node
  neighbors <- neighbors(graph, v)
  
  print(paste0("Connected Vertices to ",  Vname , ":\n"))
  
  # print(neighbors)
  
  connected_vertices <- V(graph)$preferred_name[c(neighbors)]
  
  print(connected_vertices <- unique(sort(connected_vertices)))
  
  # paths <- all_shortest_paths(graph, from = neighbor)
  
  # secondary_paths <- data.frame(neighbor = neighbor, path = unlist(paths$res))
  
  # secondary_paths %>% as_tibble() %>% distinct(neighbor, path) 
  
  return(connected_vertices)
}

shortest_paths <- lapply(QUERY, function(x) {all_neighbors_shortest_paths(graph = stringg, v = x)})

shortest_paths <- unlist(shortest_paths)

any(shortest_paths %in% QUERY) 

shortest_paths <- shortest_paths[!shortest_paths %in% QUERY]

# SCheck
any(shortest_paths %in% QUERY) # Must be false

shortest_paths <- unique(shortest_paths)


shortest_paths <- shortest_paths[!shortest_paths %in% common_vertices]


