
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

QUERYDB <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST == "CONTRAST_C") %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "110 hpf", "24 hpf"))


QUERYDB <- read_rds(list.files(path = dir, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds", full.names = T)) %>%
  filter( padj < 0.05  & abs(log2FoldChange) > 1) %>%
  filter(CONTRAST %in% c("CONTRAST_A", "CONTRAST_B")) %>%
  mutate(HPF = ifelse(sign(log2FoldChange) == -1, "pH 7.6", "pH 8.0"))

QUERYDB %>% count(CONTRAST, HPF)

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

# Module info

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"


print(WGCNA <- read_rds(paste0(wd,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))

WGCNA <- WGCNA %>% distinct(gene_id, WGCNA_target)

# First ----

# STRING
dir <- "~/Documents/MIRNA_HALIOTIS/ENSEMBLE/calcif/"

stringg <- read_rds(file.path(dir, "protein_links_full_v12.rds"))


stringg <- stringg %>% activate("nodes") %>% 
  # mutate(betweenness = betweenness(.), degree = centrality_degree(),
    # pageRank = page_rank(.)$vector) %>%
  filter(centrality_degree() > 0)

stringg <- stringg %>% activate("nodes") %>% 
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
    pageRank = page_rank(.)$vector) 
  
nodes <- stringg %>% activate("nodes") %>% as_tibble()
edges <- stringg %>% activate("edges") %>% as_tibble()

nodes %>% 
  arrange(desc(degree)) %>%
  mutate(preferred_name = factor(preferred_name, levels = unique(preferred_name))) %>%
  ggplot(aes(x = preferred_name, y = degree)) + geom_point(shape = 13) + geom_step(group = 1) +
  theme_classic(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'none', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1))

DF <- DB %>% 
  drop_na(STRINGID) %>% 
  mutate(STRINGID = strsplit(STRINGID, ";")) %>%
  unnest(STRINGID) %>% 
  mutate(SMPID = ifelse(!is.na(SMPID), "(*) ", "")) %>% 
  mutate(STRINGID_label = paste0(SMPID, STRINGID)) %>% 
  distinct(MajorRNA, STRINGID_label, STRINGID) %>%
  dplyr::count(STRINGID_label, STRINGID) %>% 
  dplyr::rename("preferred_name" = "STRINGID", "microRNA_protein" = "n") %>% 
  right_join(nodes, by = "preferred_name") %>%
  dplyr::rename("protein_protein" = "degree") %>%
  arrange(desc(protein_protein)) %>%
  mutate(preferred_name = factor(STRINGID_label, levels = unique(STRINGID_label))) %>%
  select(preferred_name, protein_protein, microRNA_protein)

# DF %>% 
#   ggplot(aes(protein_protein, microRNA_protein)) + geom_text(aes(label = preferred_name))

DF %>%
  pivot_longer(-preferred_name, values_to = "dregree") %>%
  ggplot(aes(x = preferred_name, y = dregree, color = name, group = name)) + 
  geom_point(shape = 13) + geom_step() + 
  theme_classic(base_family = "GillSans", base_size = 12) +
  theme(legend.position = 'top', 
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10))


DB <- DB %>% mutate(COG_name = ifelse(is.na(COG_name), description, COG_name))

nodes <- DB %>% 
  right_join(QUERYDB, by = "MajorRNA") %>%
  # distinct(gene_id, HPF, COG_name) %>%
  # group_by(gene_id, COG_name) %>%
  distinct(gene_id, HPF) %>%
  group_by(gene_id) %>%
  summarise(across(HPF, .fns = paste_col)) %>%
  right_join(nodes, by = "gene_id") %>%
  left_join(WGCNA, by = "gene_id")

nodes <- data.frame(
  id = as.numeric(factor(nodes$Node)), 
  label = nodes$preferred_name,
  group = nodes$HPF,
  color = nodes$WGCNA_target, # shape
  value = nodes$degree,
  nodes)

evidence <- data.frame(edges$fusion > 0,
edges$neighborhood > 0,
edges$cooccurence > 0,
edges$experiments > 0,
edges$textmining > 0,
edges$database > 0,
edges$coexpression > 0)


edges <- edges[rowSums(evidence) > 0,]
evidence <- evidence[rowSums(evidence) > 0,]

edges <- data.frame(
  value = rowSums(evidence),
  # label = paste0(rowSums(evidence), " evidences"),
  edges
)

# edges$from <- nodes[edges$from,]$Node
# edges$to <- nodes[edges$to,]$Node
library(visNetwork)

visNetwork::visNetwork(
  nodes = nodes,
  edges = edges, width = "100%") %>%
  visInteraction(navigationButtons = T) %>%
  visOptions(highlightNearest = T, nodesIdSelection = T, selectedBy = "COG_name") %>%
  visLegend(
    useGroups = FALSE, addNodes = data.frame(label = "Protein", shape = "circle"), 
    addEdges = data.frame(label = "Predicted association", color = "black"))
  # visEdges(arrows = 'from')

# visNetwork(nodes, edges, width = "100%") %>% 
#   visEdges(arrows = "to") %>% 
#   visHierarchicalLayout(direction = "DU") 

# Second ----

# Third (hard to try)
# split miRs to viz
qnodes <- DB %>% right_join(QUERYDB) %>% 
  distinct(gene_id, MajorRNA) %>%
  count(gene_id, sort = T) %>% sample_n(2)

q <- qnodes%>% pull(gene_id)

stringg %>% activate("nodes") %>% filter(gene_id %in% q)
  mutate(Col = ifelse(preferred_name %in% which_targets, "Target under OA", ""),
    label = ifelse(preferred_name %in% which_targets, preferred_name, ""))



Net <- DB %>% right_join(QUERYDB) %>% filter(gene_id %in% q)

# str(qnodes <- QUERYDB %>% pull(MajorRNA))
# Net <- DB %>% filter(MajorRNA %in% qnodes)

Net <- Net %>% distinct(MajorRNA, gene_id) %>%
  dplyr::rename("from" = "MajorRNA", "to" = "gene_id")


NodeInfo <- DB %>% 
  distinct(gene_id, COG_category) %>% 
  drop_na(COG_category) %>%
  group_by(gene_id) %>% 
  sample_n(1) %>%
  rename("Node" = "gene_id")

NodeInfo2 <- QUERYDB %>% 
  distinct(MajorRNA, HPF, WGCNA, MirGeneDB_ID) %>% 
  rename("Node" = "MajorRNA")

# Make it hierarchical

# create a data frame giving the hierarchical structure of your individuals

root <- data.frame(from = "root", to = unique(Net$from))

nrow(hierarchy <- rbind(root, Net))

hierarchy <- hierarchy %>% as_tibble()

hierarchy <- hierarchy %>% 
  left_join(NodeInfo2, by = c("to" = "Node"))

# view(hierarchy)

# hierarchy %>% mutate()
# create a vertices data.frame. One line per object of our hierarchy

Nodes <- data.frame(Node = unique(c(hierarchy$from, hierarchy$to))) 

Nodes%>% 
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

# Use class inheritance for layout but plot class imports as bundles
edge_bundle_plot <- ggraph(graph, 'dendrogram', circular = T) +
  theme_void()

edge_bundle_plot$data["angle"] <- sapply(1:nrow(edge_bundle_plot$data), 
  function (i) {
    t_ratio <- min(max(edge_bundle_plot$data$y[i] / edge_bundle_plot$data$x[i], -10000),10000)
    atan(t_ratio) * 180 / pi
  })

edge_bundle_plot$data["hjust"] <- sapply(edge_bundle_plot$data$x, 
  function (x) {
    ifelse(x < 0, 1, 0)
  })

edge_bundle_plot <- edge_bundle_plot +
  geom_conn_bundle(aes(colour = after_stat(index)),
    data = get_con(from, to),
    edge_alpha = 0.25,
    width=0.9, 
    tension=1
  ) +
  coord_fixed()

# Create Nodes
edge_bundle_plot +
  geom_node_point(aes(filter = leaf, x = x, y=y), alpha=0.8) +
  geom_node_text(aes(filter = leaf, x = x*1.05, y=y*1.05, 
    label=Node, # colour=region, alpha=ranking,
    angle=angle, hjust=hjust), size=2.7) 



# Basic graph
# https://github.com/thomasp85/ggraph/issues/272
# https://lemuelkumarga.com/viewer/?key=wildlife-trade
# https://r-graph-gallery.com/310-custom-hierarchical-edge-bundling.html



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


ggraph(graph, 'matrix', sort.by = node_rank_spectral()) + 
  geom_edge_point(mirror = TRUE) + 
  coord_fixed()

ggraph(graph, 'stress') + 
  geom_node_voronoi(aes(fill = WGCNA), max.radius = 0.05) + 
  # geom_node_point() + 
  geom_edge_link() + 
  coord_fixed()


# Find relation between target net and stringnet -------
# Subst only to degs

m1 <- DB %>%
  distinct(MajorRNA, STRINGID) %>%
  drop_na(STRINGID) %>%
  # filter(STRINGID %in% hc2$labels) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = MajorRNA, values_from = n, values_fill = 0) %>%
  data.frame() 

rownames(m1) <- m1$STRINGID

m1$STRINGID <- NULL

g <- stringg

nodes <- g %>% activate("nodes") %>% as_tibble()
edges <- g %>% activate("edges") %>% as_tibble()

edges$from <- nodes[edges$from,]$preferred_name
edges$to <- nodes[edges$to,]$preferred_name

m2 <- edges %>% 
  drop_na() %>%
  select(from, to) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = from, values_from = n, values_fill = 0) %>%
  data.frame() 

rownames(m2) <- m2$to

m2$to <- NULL


# filter matrix

dim(m1)
dim(m2)

sum(keep <- rownames(m1) %in% rownames(m2))

m1 <- m1[keep,]

sum(keep <- rownames(m2) %in% rownames(m1))

m2 <- m2[keep,]


hc1 <- stats::hclust(dist(m1, method = "binary"), method="ward.D2")

hc2 <- stats::hclust(dist(m2, method = "binary"), method="ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

# Align and plot two dendrograms side by side
dendlist("targetnet" = dend1, "stringnet" = dend2) %>%
  untangle(method = "random") %>% # Find the best alignment layout
  tanglegram(
    sort = T,
    highlight_distinct_edges = F,
    highlight_branches_col = F,
    highlight_branches_lwd = F,
    common_subtrees_color_lines = F,
    common_subtrees_color_branches = T) # %>% 
  # plot(sub = paste("entanglement =", round(entanglement(.), 2)))


cor_bakers_gamma(dend1, dend2)


dendlist(dend1, dend2) %>%
  untangle(method = "random", R = 20) %>%
  entanglement() # Alignment quality

# Create a list to hold dendrograms
dend_list <- dendlist(dend1, dend2)

cors <- cor.dendlist(dend_list, method = "common_nodes")
# Print correlation matrix
round(cors, 2)
