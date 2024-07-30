
library(tidyverse)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

QUERYDB <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "24 hpf", "110 hpf"))

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

str(qnodes <- QUERYDB %>% pull(MajorRNA))

Net <- DB %>% filter(MajorRNA %in% qnodes)

Net <- Net %>% distinct(MajorRNA, gene_id) %>%
  dplyr::rename("from" = "MajorRNA", "to" = "gene_id")


NodeInfo <- DB %>% 
  distinct(gene_id, COG_category) %>% 
  drop_na(COG_category) %>%
  group_by(gene_id) %>% 
  sample_n(1) %>%
  rename("Node" = "gene_id")

NodeInfo2 <- QUERYDB %>% 
  distinct(MajorRNA, HPF, WGCNA) %>% 
  rename("Node" = "MajorRNA")

# Make it hierarchical

# create a data frame giving the hierarchical structure of your individuals

root <- data.frame(from = "origin", to = unique(Net$from))

hierarchy <- rbind(root, Net)

hierarchy <- hierarchy %>% as_tibble()

hierarchy <- hierarchy %>% left_join(NodeInfo2, by = c("from" = "Node"))

# hierarchy %>% mutate()
# create a vertices data.frame. One line per object of our hierarchy

Nodes <- data.frame(Node = unique(c(hierarchy$from, hierarchy$to))) %>% 
  mutate_if(is.character, as.factor) 

Nodes <- Nodes %>% left_join(NodeInfo2) %>% left_join(NodeInfo)

# Net %>% count(to) %>% Nodes

# Let's add a column with the group of each name. It will be useful later to color points
# vertices$group  <-  hierarchy$from[ match( vertices$Node, hierarchy$to ) ]

# graph <- graph_from_data_frame( hierarchy, vertices=Nodes )

graph <- tidygraph::tbl_graph(nodes = Nodes, edges = hierarchy, directed = T)


from <- match(Net$from, Nodes$Node)
to <- match(Net$to, Nodes$Node)


# graph %>% activate("nodes") %>% filter(!is.na(COG_category))

graph <- graph %>% activate("nodes") %>%
  mutate(betweenness = betweenness(.), degree = centrality_degree())

# Basic graph
# https://github.com/thomasp85/ggraph/issues/272
# https://lemuelkumarga.com/viewer/?key=wildlife-trade
# https://r-graph-gallery.com/310-custom-hierarchical-edge-bundling.html


# Use class inheritance for layout but plot class imports as bundles
ggraph(graph, 'dendrogram', circular = TRUE) +
  # geom_conn_bundle(aes(colour = HPF, group = HPF),
  geom_conn_bundle(aes(colour = after_stat(index)),
    data = get_con(from, to),
    edge_alpha = 0.25,
    width=0.9, 
    tension=1
  ) +
  geom_node_point(aes(filter = leaf, colour = COG_category)) +
  scale_edge_colour_distiller('', direction = 1, guide = 'edge_coloursteps') +
  coord_fixed() +
  ggforce::theme_no_axes()
  # facet_nodes(~ HPF)


# Continue w/ other netviz


# Binary model matrix
# one category
# x <- factor(Net$MajorRNA)
# m <- model.matrix(~ x )
# attr(m, "ATT") <- NULL

# heatmap(m)

# multiple category
# x <- factor(Net$MajorRNA)
# y <- factor(Net$gene_id)
# m <- model.matrix(~ x + y -1)
# diag(m) <- 0

# g <- igraph::graph.adjacency(m)




Nodes <- data.frame(Node = unique(c(Net$from,Net$to))) %>% 
  mutate_if(is.character, as.factor) 

Nodes <- Nodes %>% left_join(NodeInfo) %>% left_join(NodeInfo2)

# Edges <- data.frame(from = net$gene_id, to = net$MajorRNA)
graph <- tidygraph::tbl_graph(nodes = Nodes, edges = Net, directed = T)

# graph <- graph_from_data_frame(Nodes, vertices=Net )


library(igraph)
library(tidygraph)
library(ggraph)

graph <- graph %>% activate("nodes") %>%
  mutate(betweenness = betweenness(.), degree = centrality_degree())

# igraph::E(graph)$weight$membership

# igraph::communities(igraph::cluster_infomap(graph))

graph <- graph %>% activate("nodes") %>%
  mutate(
    # infomap = igraph::cluster_infomap(.),
    pageRank = page_rank(.)$vector)

layout <- create_layout(graph, layout = 'igraph', algorithm = 'kk')


ggraph(layout) +
  geom_edge_arc(
    strength = 0.1,
    start_cap = circle(1.5, 'mm'),
    end_cap = circle(1.5, 'mm'),
    arrow = arrow(
      angle = 90,
      length = unit(0.05, "inches"),
      ends = "first",
      type = "open")) +
  geom_node_point(aes(color = WGCNA)) + 
  # geom_node_text(aes(label = Node), repel = T, size = 1.5) +
  coord_fixed() +
  theme_graph(base_family = "GillSans")


ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_node_point(aes(color = WGCNA)) 


ggraph(graph, 'matrix', sort.by = node_rank_spectral()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed()
