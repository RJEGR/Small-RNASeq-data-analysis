
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

RES.p <- RES %>% 
  filter(CONTRAST %in% c("CONTRAST_A","CONTRAST_B")) %>%
  # mutate(CONTRAST = dplyr::recode_factor(CONTRAST, !!!recode_to)) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

QUERIES <- RES.p %>% 
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
  mutate(OA_MIR_degree) %>%
  mutate(OA_MIR_degree = ifelse(is.na(OA_MIR_degree), 0, OA_MIR_degree))
  # dplyr::rename("preferred_name" = "STRINGID")


NodeDF <- DF2 %>% right_join(DF1 , by = "gene_id") 

NodeDF <- DB %>% 
  # mutate(COG_name = ifelse(is.na(COG_name), "Unknown", STRINGID)) %>%
  mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
  distinct(gene_id, COG_name) %>%
  mutate(COG_name = factor(COG_name, levels = sort(decreasing = F, unique(COG_name)))) %>%
  right_join(NodeDF, by = "gene_id")

print(BT <- read_tsv(file.path(dir, "Supplementary_tables - Biological_themes.tsv"), col_names = T))

BT <- BT %>% select(-microRNA_degree, -STRINGID)

NodeDF <- NodeDF %>%
  right_join(BT, by = "gene_id")

view(NodeDF)
# Plot something intersting
# as protein

DataViz <- NodeDF %>%
  mutate(STRINGID = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>% 
  filter(!COG_name %in% "NA, NA") %>%
  mutate(Frac = OA_MIR_degree) %>%
  filter(Frac > 0) %>%
  group_by(Contrast) %>%
  arrange(desc(Frac), .by_group = T) %>% 
  mutate(Label = STRINGID, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
        levels = rev(paste(Label, row_number, sep = "__")))) %>%
  mutate(Contrast = factor(Contrast, levels = c("24 hpf", "110 hpf")))

DataViz %>%
  # mutate(Label = factor(Label, levels = unique(Label))) %>%
  ggplot(aes(y = Label, x = OA_MIR_degree, color = Contrast)) +
  # facet_grid( COG_name ~ ., space = "free_y", scales = "free_y") +
  ggforce::facet_col(~ COG_name, space = "free", scales = "free_y") +
  # facet_wrap(~ COG_name, ncol = 1, scales = "free_y") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  scale_color_manual(values = c("black", "gray")) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "gray",hjust = 0),
    legend.position = 'top',
    panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) -> P




dat_text <- DataViz %>% 
  # select(Label) %>%
  mutate(Text = gsub("__.+$", "", Label), 
        Text = factor(Text))
  
P <- P +
  geom_text(data = dat_text, 
    aes(label=Text), hjust= -0.25, vjust = 0, size = 1.5, family = "GillSans", color = "gray20") +
  scale_x_continuous("MicroRNA degree", limits = c(0,10), breaks = seq(0,10, by = 2))

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

ggsave(P, filename = 'Mirdegree_COGS.png', path = dir, width = 4, height = 10, device = png, dpi = 500)

# DataViz %>%
#   ggplot(aes(x = ))


DataVizsum <- DataViz %>%
  mutate(COG_name = as.character(COG_name)) %>%
  group_by(Contrast, COG_name) %>%
  summarise(OA_MIR_degree = sum(OA_MIR_degree), MIR_degree = sum(MIR_degree)) %>%
  group_by(Contrast) %>%
  arrange(desc(OA_MIR_degree), .by_group = T) %>% 
  mutate(Label = as.character(COG_name), row_number = row_number(Label)) %>% 
  # mutate(Label = paste(Label, row_number, sep = "__")) %>%
  mutate(Label = factor(paste(Label, row_number, sep = "__"), levels = rev(paste(Label, row_number, sep = "__")))) 


DataVizsum %>%
  ggplot(aes(y = Label, x = OA_MIR_degree)) +
  # facet_grid( COG_name ~ ., space = "free_y", scales = "free_y") +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  # facet_wrap(~ COG_name, ncol = 1, scales = "free_y") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("black", "black")) +
  # scale_x_continuous("MicroRNA degree", limits = c(0,20), breaks = seq(0,20, by = 5)) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    legend.position = 'top',
    # panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) -> P



P

ggsave(P, filename = 'Mirdegree_sum_COGS.png', path = dir, width = 4.5, height = 4, device = png, dpi = 500)

# fold-increase decrease 

RESviz <- RES.p %>% distinct(MajorRNA, MirGeneDB_ID, log2FoldChange) %>%
  # mutate(Contrast = ifelse(MirGeneDB_ID == "MIR-133-3p" & log2FoldChange > 0, NA, "")) %>% drop_na(Contrast) %>%
  filter(log2FoldChange < 0) %>%
  mutate(Contrast = ifelse(MirGeneDB_ID %in% names(QUERIES[1:3]), "24 hpf", "110 hpf")) %>%
  mutate(Contrast = factor(Contrast, levels = c("24 hpf", "110 hpf")))

# Lev110hpf <- DataVizsum %>% ungroup %>% filter(Contrast %in% "110 hpf") %>% distinct(COG_name) %>% pull()


DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
  filter(!COG_name %in% "NA, NA") %>%
  distinct(MajorRNA, COG_name) %>%
  # left_join(DataVizsum %>% ungroup() %>% distinct(COG_name, row_number))
  right_join(RESviz, by = "MajorRNA")

DF1 <- DataVizsum %>% 
  mutate(names_to = "OA_MIR_degree") %>%
  dplyr::rename("value_to" = "OA_MIR_degree") %>%
  select(Contrast, Label, names_to, value_to)

DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
  filter(!COG_name %in% "NA, NA") %>%
  distinct(MajorRNA, COG_name) %>%
  # left_join(DataVizsum %>% ungroup() %>% distinct(COG_name, row_number))
  right_join(RESviz, by = "MajorRNA") %>%
  mutate(Label = paste(COG_name, Contrast, sep = "__")) %>%
  mutate(names_to = "log2FoldChange") %>%
  dplyr::rename("value_to" = "log2FoldChange") %>%
  select(Contrast, Label, names_to, value_to) %>%
  rbind(DF1) %>%
  ggplot(aes(y = Label, x = value_to)) +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  # scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  # geom_jitter(alpha = 0.5) +
  # stat_boxplot(geom ='errorbar', width = 0.3, linetype="dashed") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", alpha = 0.5) +
  facet_grid( Contrast ~ names_to, space = "free", scales = "free") 


DB %>% 
  filter(MajorRNA %in% QUERIES) %>%
  mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
  filter(!COG_name %in% "NA, NA") %>%
  distinct(MajorRNA, COG_name) %>%
  # left_join(DataVizsum %>% ungroup() %>% distinct(COG_name, row_number))
  right_join(RESviz, by = "MajorRNA") %>%
  mutate(Label = paste(COG_name, Contrast, sep = "__")) %>%
  # filter(COG_name %in% as.character(levelshpf)) %>%
  mutate(Label = factor(Label, levels = rev(levelshpf))) %>%
  ggplot(aes(y = Label, x = log2FoldChange)) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) +
  # geom_jitter(alpha = 0.5) +
  stat_boxplot(geom ='errorbar', width = 0.3, linetype="dashed") +
  geom_point(shape = 21, size = 1, stroke = 1.2, fill = "white", alpha = 0.5) +
  # geom_segment(aes(xend = 0, yend = COG_name), linewidth = 0.5, color="#E7DFD5") +
  # geom_boxplot(width = 0.3, outlier.alpha = 0, linetype="dashed") +
  ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # scale_color_manual(values = c("black", "black")) +
  # scale_x_continuous("MicroRNA degree", limits = c(0,20), breaks = seq(0,20, by = 5)) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 0),
    legend.position = 'top',
    # panel.border = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    panel.grid.minor = element_blank(),
    # axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title = element_blank(),
    panel.grid.major = element_blank()) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) -> P2

# P2 +  geom_text(data = RESviz, 
#   aes(label=MirGeneDB_ID), hjust= -0.25, vjust = 0, size = 1.5, family = "GillSans", color = "gray20") 
# 

# DB %>% 
#   filter(MajorRNA %in% QUERIES) %>%
#   mutate(COG_name = paste0(COG_category, ", ",COG_name)) %>%
#   filter(!COG_name %in% "NA, NA") %>%
#   distinct(MajorRNA, COG_name) %>%
#   # left_join(DataVizsum %>% ungroup() %>% distinct(COG_name, row_number))
#   right_join(RESviz, by = "MajorRNA") %>%
#   mutate(COG_name = paste(COG_name, Contrast, sep = "__")) %>%
#   filter(COG_name %in% as.character(levelshpf)) %>%
#   mutate(COG_name = factor(COG_name, levels = rev(as.character(levelshpf)))) %>%
#   ggplot(aes(y = COG_name, x = log2FoldChange)) +
#   # ggforce::facet_col(~ Contrast, space = "free", scales = "free_y") +
#   ggridges::geom_density_ridges_gradient(
#   jittered_points = T,
#   position = ggridges::position_points_jitter(width = 0.05, height = 0),
#   point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) 
#   # scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
# 


library(patchwork)

ps <- P + plot_spacer() + P2 + plot_layout(widths = c(5, -0.1, 4))

ps

ggsave(ps, filename = 'WGCNA2NOGS.png', path = dir, width = 6, height = 4, device = png, dpi = 300)


# Visualice network

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


