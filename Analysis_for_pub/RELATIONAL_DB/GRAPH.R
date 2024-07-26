
rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))


qnodes <- QUERYDB %>% filter(HPF == "24 hpf") %>% pull(MajorRNA)

qnodes <- QUERYDB %>% group_by(HPF) %>% sample_n(10) %>% pull(MajorRNA)

Net <- DB %>% filter(MajorRNA %in% qnodes)

NodeInfo <- Net %>% distinct(gene_id, COG_category) %>% rename("Node" = "gene_id")

NodeInfo2 <- QUERYDB %>% distinct(MajorRNA, HPF) %>% 
  rename("Node" = "MajorRNA")

# Binary model matrix
# one category
x <- factor(Net$MajorRNA)
m <- model.matrix(~ x )
attr(m, "ATT") <- NULL

heatmap(m)

# multiple category
# x <- factor(Net$MajorRNA)
# y <- factor(Net$gene_id)
# m <- model.matrix(~ x + y -1)
diag(m) <- 0
g <- igraph::graph.adjacency(m)

Net <- Net %>% distinct(gene_id, MajorRNA)

Nodes <- data.frame(Node = unique(c(Net$MajorRNA,Net$gene_id))) %>% 
  mutate_if(is.character, as.factor) 

Nodes <- Nodes %>% left_join(NodeInfo) %>% left_join(NodeInfo2)

# Edges <- data.frame(from = net$gene_id, to = net$MajorRNA)
graph <- tidygraph::tbl_graph(nodes = Nodes, edges = Net, directed = T)

library(igraph)
library(tidygraph)
library(ggraph)

graph <- graph %>% activate("nodes") %>%
  mutate(betweenness = betweenness(.), degree = centrality_degree())

# igraph::E(graph)$weight$membership

igraph::communities(igraph::cluster_infomap(graph))

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
  geom_node_point(aes(color = HPF)) + 
  # geom_node_text(aes(label = Node), repel = T, size = 1.5) +
  coord_fixed() +
  theme_graph(base_family = "GillSans")


ggraph(graph, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(
    strength = 0.1,
    start_cap = circle(1.5, 'mm'),
    end_cap = circle(1.5, 'mm'),
    arrow = arrow(
      angle = 90,
      length = unit(0.05, "inches"),
      ends = "first",
      type = "open")) +
  geom_node_point(aes(color = HPF)) 
