
# Dominant regulatory roles of network hubs. The tar-
# Dominant regulatory roles of network hubs. The tar- geting of many genes by each miRNA, the targeting of individual genes by multiple miRNAs and the down- stream effects that result from the miRNA-mediated regulation of transcription factors lead to highly com- plex networks of miRNAs and their target genes34 (FIG. 1). The nodes of these networks (which can be individual miRNAs or mRNAs, including transcription factors), are typically connected to many other nodes in the complex regulatory webs. 

# It is useful to identify nodes with atypi- cally high numbers of connections (‘hubs’) because these represent sites of signalling convergence with potentially large explanatory power for network behaviour 

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}



dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

gene_read_exp_df <- read_tsv(file.path(dir, "Transcriptome_read_expression.tsv"))


gene_read_exp_df %>% 
  group_by(dpf,gene_id) %>%
  summarise(Read_expression = sum(Read_expression)) %>%
  ggplot(aes(y = gene_id, x = as.factor(dpf), fill = log10(Read_expression))) + geom_tile()


gene_read_exp_df <- gene_read_exp_df %>% 
  group_by(gene_id) %>%
  summarise(Read_expression = sum(Read_expression)) 

# Tag proteins important for acidification
# PREPARING NODE SHAPE 

QUERYDB <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "pH 7.6", "pH 8.0"))

QUERYDB %>% count(CONTRAST, HPF)

QUERYMIRS <- QUERYDB %>% distinct(MajorRNA) %>% pull()

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")) %>% 
    mutate(STRINGID = strsplit(STRINGID, ";")) %>%
    unnest(STRINGID) %>%
    mutate(Contrast = strsplit(Contrast, ";")) %>%
    unnest(Contrast))


paste_col <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

# Using only up-expressed by acidification

recode_to <- c(`CONTRAST_B_Low`= "110 hpf under 7.6", 
  `CONTRAST_A_Low` = "24 hpf under 7.6", 
  `CONTRAST_A_Low;CONTRAST_B_Low` = "Both under pH 7.6"
  )

#Include

QTARGETSDB <- DB %>% filter(MajorRNA %in% QUERYMIRS) %>% 
  filter(grepl("CONTRAST_B_Low|CONTRAST_A_Low", Contrast)) %>%
  distinct(Contrast,STRINGID) %>%
  drop_na() %>% 
  group_by(STRINGID) %>%
  summarise(across(Contrast, .fns = paste_col), n = n()) %>% 
  # filter(n == 1) %>%
  dplyr::mutate(HPF = dplyr::recode_factor(Contrast, !!!recode_to)) %>%
  # dplyr::mutate(PH = "pH 7.6") %>%
  dplyr::select(-n) %>%
  ungroup()

QTARGETSDB %>% count(HPF)

# QTARGETS <- DB %>% filter(MajorRNA %in% QUERYMIRS) %>% distinct(STRINGID) %>% drop_na() %>% pull()

print(BT <- read_tsv(file.path(dir, "Supplementary_tables - Biological_themes.tsv"), col_names = T))

BT <- BT %>% arrange(desc(microRNA_degree))

# Processing biological themes Database

BT <- BT %>% 
  drop_na(Bilogical_pathway) %>% 
  mutate(Bilogical_pathway = strsplit(Bilogical_pathway, ",")) %>%
  unnest(Bilogical_pathway) %>% 
  # mutate(SMPID = ifelse(!is.na(SMPID), "(*) ", "")) %>% 
  # mutate(STRINGID_label = paste0(SMPID, STRINGID)) %>% 
  # distinct(MajorRNA, STRINGID_label, STRINGID) %>%
  dplyr::count(gene_id, sort = T) %>% 
  dplyr::rename("N_pathways" = "n") %>% right_join(BT)

BT <- BT %>% mutate(Bilogical_pathway = ifelse(N_pathways > 1, "Diverse_pathway", Bilogical_pathway))

BT <- BT %>% left_join(QTARGETSDB) %>% 
  mutate(HPF = ifelse(is.na(HPF), "Basal", as.character(HPF)))

BT %>% count(HPF)

dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

print(stringg <- read_rds(file.path(dir, "protein_links_full_v12.rds")))

# Processing edge evidence

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

# BT <- BT %>% distinct(gene_id, microRNA_degree, Bilogical_pathway, HPF)
# Bind info to nodes

stringg <- stringg %>% activate("nodes") %>% left_join(gene_read_exp_df, by = "gene_id") 

stringg <- stringg %>% activate("nodes") %>%  left_join(BT, by = "gene_id") 

# stringg <- stringg %>% activate("nodes") %>% as_tibble()

# Get edges connected to the specified node ----

# QUERY <- c("MAP11") # RNF166

QUERY <- QTARGETSDB %>% distinct(STRINGID) %>% pull()


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

# common_vertices <- common_vertices[!common_vertices %in% QUERY]

# If decide walk the last path -----

# Initialize a list to store secondary paths
secondary_paths <- list()

# Loop through each neighbor and find paths to other nodes
for (neighbor in common_vertices) {
  
  cat("Conected: ", neighbor, "\n")
  
  # Find all shortest paths from the neighbor to all other nodes
  paths <- all_shortest_paths(stringg, from = neighbor)
  
  # secondary_paths[[as.character(neighbor)]] <- paths$res
  
  secondary_paths[[as.character(neighbor)]] <- data.frame(neighbor = neighbor, path = unlist(paths$res))
  
}

# Create a data frame
secondary_paths <- do.call(rbind, secondary_paths) %>% as_tibble() %>% distinct(neighbor, path) 

secondary_paths %>% group_by(neighbor) %>% dplyr::count()

secondary_paths <- secondary_paths %>% group_by(neighbor) %>% slice_head(n = 10)

str(connected_vertices <- unique(sort(secondary_paths$path)))

# MAP <- secondary_paths
# 
# MAP <- data.frame(Name = rep(names(MAP),
#   lapply(MAP, length)),
#   GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble()


# Then ----
# Filter these connecte vertices 

# print(connected_vertices <- V(stringg)$preferred_name[c(node_id, connected_vertices)])

# connected_vertices <- unique(sort(connected_vertices))

as_adjacency_matrix(stringg, "upper")

layout <- stringg %>% activate("nodes") %>% 
  filter(degree > 0) %>%
  # filter(preferred_name %in% QUERY[1]) %>%
  mutate(colNode = "Basal") %>%
  mutate(colNode = ifelse(preferred_name %in% QUERY, "Target", colNode)) %>%
  mutate(colNode = ifelse(preferred_name %in% shortest_paths, "short_paths_target", colNode)) %>%
  mutate(colNode = ifelse(preferred_name %in% common_vertices, "Common_vertices", colNode)) %>%
  activate("edges") %>% 
  mutate(Evidence_size = ifelse(Evidence_size > 0.15, Evidence_size, NA)) %>%
  # filter(Evidence_size > 0.2) %>%
  create_layout(layout = 'auto')
  # create_layout(layout = 'linear')


ggraph(layout) +
  geom_edge_link(aes(alpha = Evidence_size)) +
  geom_node_point(aes(color = Bilogical_pathway, size = microRNA_degree+10, shape = colNode)) +
  geom_node_text(aes(label = preferred_name), repel = TRUE, family = "GillSans") +
  ggthemes::scale_color_calc() +
  # facet_nodes(~ Bilogical_pathway, scales = "free") +
  theme_bw(base_family = "GillSans") 
  

ggraph(stringg, layout = "matrix", sort.by = node_rank_traveller()) + 
  geom_edge_point() +
  geom_node_text(aes(label = preferred_name), repel = TRUE, family = "GillSans") 


# Deal with spagetti plot

# Calculate various network properties, adding them as attributes
# to each node/vertex

# stringg <- stringg %>% activate("nodes") %>% 
#   mutate(betweenness = betweenness(.), degree = centrality_degree(),
#     pageRank = page_rank(.)$vector) 

graph <- stringg

# take a while calculating membership

V(graph)$comm <- membership(optimal.community(graph))

# V(graph)$closeness <- centralization.closeness(graph)$res
# V(graph)$betweenness <- centralization.betweenness(graph)$res
# V(graph)$eigen <- centralization.evcent(graph)$vector


# Re-generate dataframes for both nodes and edges, now containing
# calculated network attributes

node_list <- get.data.frame(graph, what = "vertices")


names(node_list)
names(node_list)[3] <- "name"

# Determine a community for each edge. If two nodes belong to the
# same community, label the edge with that community. If not,
# the edge community value is 'NA'

edge_list <- get.data.frame(graph, what = "edges") # %>%
  # inner_join(node_list %>% select(name, comm), by = c("from" = "name")) %>%
  # inner_join(node_list %>% select(name, comm), by = c("to" = "name")) %>%
  # mutate(group = ifelse(comm.x == comm.y, comm.x, NA) %>% factor())


edge_list$from <- node_list[edge_list$from,]$name
edge_list$to <- node_list[edge_list$to,]$name

# Create a character vector containing every node name
all_nodes <- unique(sort(node_list$name))

# Adjust the 'to' and 'from' factor levels so they are equal
# to this complete list of node names
plot_data <- edge_list %>% mutate(
  to = factor(to, levels = all_nodes),
  from = factor(from, levels = all_nodes))


name_order <- (node_list %>% arrange(comm))$name

# Reorder edge_list "from" and "to" factor levels based on
# this new name_order
plot_data <- edge_list %>% mutate(
  to = factor(to, levels = name_order),
  from = factor(from, levels = name_order))


# plot
query_nodes <- c("")
query_nodes <- paste(query_nodes, collapse = "|")
links_ <- c("from", "to")
plot_data %>% filter_at(vars(all_of(links_)), any_vars(grepl(query_nodes, .)))

plot_data %>%
  filter_all()

# Create the adjacency matrix plot
ggplot(plot_data, aes(x = from, y = to)) +
  geom_raster() +
  theme_bw() +
  # Because we need the x and y axis to display every node,
  # not just the nodes that have connections to each other,
  # make sure that ggplot does not drop unused factor levels
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme(
    # Rotate the x-axis lables so they are legible
    axis.text.x = element_text(angle = 270, hjust = 0),
    # Force the plot into a square aspect ratio
    aspect.ratio = 1,
    # Hide the legend (optional)
    legend.position = "none")
