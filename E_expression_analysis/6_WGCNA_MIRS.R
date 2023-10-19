# WGCNA FOR ONLY MIRS
# 

# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"
# read_tsv(paste0(wd, "DESEQ2SRNATARGET.tsv"))

library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

head(COUNT <- read_rds(paste0(path, "COUNT.rds"))) # ONLY MIRS

DB <- read_tsv(paste0(path, "DESEQ_RES.tsv")) %>% distinct(Name, Family)

# MIRGENEDB <- read_tsv(paste0(path, "SRNA2MIRGENEDB.tsv")) %>% distinct(Name, Family, arm)
# MIRGENEDB <- MIRGENEDB %>% mutate(Fam = Family) %>% tidyr::unite("MIR", c("Fam", "arm"), sep = "-") 
# DB <- DB %>% left_join(MIRGENEDB)


# COLLAPSE BY FAM

DB %>% dplyr::count(Name, sort = T)

COUNT <- COUNT %>% as_tibble(rownames = "Name") %>%
  right_join(DB %>% distinct(Name, Family), by = "Name") %>%
  group_by(Family) %>%
  summarise_at(vars(colnames(COUNT)), sum) 

which_cols <- COUNT %>% dplyr::select_if(is.double) %>% colnames()

datExpr <- COUNT %>% dplyr::select(all_of(which_cols))

datExpr <- as(datExpr, "matrix")

rownames(datExpr) <- COUNT$Family

# THEN,

# datExpr <- log2(datExpr+1)

datExpr <- t(datExpr) # log2(count+1) # 

str(datExpr)

cat("\n:::::\n")

gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

max_power <- 20

powers = c(c(1:10), seq(from = 10, to = max_power, by=1)) 

#powers = unique(powers)

allowWGCNAThreads()

cor_method =  "cor" # by default WGCNA::cor(method =  'pearson') is used, "bicor"

corOptionsList = list(use ='p') # maxPOutliers = 0.05, blocksize = 20000

sft = pickSoftThreshold(datExpr, 
  powerVector = powers, 
  corFnc = cor_method,
  corOptions = corOptionsList,
  verbose = 5, 
  networkType = "unsigned")

# sft <- read_rds(paste0(path, 'SoftThreshold_',cor_method, '.rds'))

soft_values <- abs(sign(sft$fitIndices[,3])*sft$fitIndices[,2])

soft_values <- round(soft_values, digits = 2)

hist(soft_values)

power_pct <- quantile(soft_values, probs = 0.95)

softPower <- sft$fitIndices[,1][which(soft_values >= power_pct)]

meanK <- sft$fitIndices[softPower,5]

hist(sft$fitIndices[,5])

softPower <- min(softPower)

cat("\nsoftPower value", softPower, '\n')


title1 = 'Scale Free Topology Model Fit,signed R^2'
title2 = 'Mean Connectivity'

caption = paste0("Lowest power for which the scale free topology index reaches the ", power_pct*100, " %")

sft$fitIndices %>% 
  mutate(scale = -sign(slope)*SFT.R.sq) %>%
  select(Power, mean.k., scale) %>% pivot_longer(-Power) %>%
  mutate(name = ifelse(name %in% 'scale', title1, title2)) %>%
  ggplot(aes(y = Power, x = value)) +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  facet_grid(~name, scales = 'free_x', switch = "x") +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave

# Build a adjacency "correlation" matrix

enableWGCNAThreads()

# specify network type

adjacency <- adjacency(datExpr, 
  power = softPower, 
  type = "unsigned")


TOM <- TOMsimilarity(adjacency, TOMType = "signed") # specify network type


dissTOM = 1 - TOM

heatmap(TOM)

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", 
  main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module

minClusterSize <- 5

dynamicMods <- cutreeDynamic(dendro= geneTree, 
  distM = dissTOM,
  method = "hybrid",
  deepSplit = 2, 
  cutHeight = 0.8,
  pamRespectsDendro = FALSE,
  minClusterSize = minClusterSize)

dynamicColors = labels2colors(dynamicMods)
names(dynamicColors) <- colnames(datExpr)
MEList = moduleEigengenes(datExpr, 
  colors= dynamicColors,
  softPower = softPower)

MEs = MEList$eigengenes

MEDiss= 1 - cor(MEs)

METree = flashClust(as.dist(MEDiss), method= "average")

# Set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0

MEDissThres = 0.0

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

#plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic\n(cutHeight=0.2)"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))

print(nm)

write_rds(as_tibble(moduleColors, rownames = "Family") %>% rename("value" = "WGCNA"), 
  file = paste0(path, "WGCNA_MIRS.rds"))

# Plot TOM
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
TOMplot(plotTOM, dendro = geneTree, Colors = moduleColors)

#

reads <- colSums(datExpr)

Total <- sum(reads)

# Sanity check:

identical(names(colSums(datExpr)), names(moduleColors))

data.frame(reads, moduleColors) %>% 
  as_tibble(rownames = "Name") %>% 
  group_by(moduleColors) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'moduleColors') -> stats

# Prep binary datTraits

datTraits <- gsub(".clean.newid.subset", "", rownames(datExpr))

HR11076 <- grepl('HR11076', datTraits)

HR1108 <- grepl('HR1108', datTraits)

HR2476 <- grepl('HR2476', datTraits)

HR248 <- grepl('HR248', datTraits)

datTraits <- data.frame(HR11076, HR1108, HR2476, HR248)

datTraits <- 1*datTraits

rownames(datTraits) <- rownames(datExpr)

# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

MEs = orderMEs(MEs0)

names(MEs) <- str_replace_all(names(MEs), '^ME', '')

moduleTraitCor = cor(MEs, datTraits, use= "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datTraits))


moduleTraitCor %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'moduleTraitCor') -> df1

moduleTraitPvalue %>% as_tibble(rownames = 'module') %>% 
  pivot_longer(-module, values_to = 'corPvalueStudent') %>%
  right_join(df1) -> df1

hclust <- hclust(dist(moduleTraitCor), "complete")

df1 %>%
  mutate(star = ifelse(corPvalueStudent <.001, "***", 
    ifelse(corPvalueStudent <.01, "**",
      ifelse(corPvalueStudent <.05, "*", "")))) -> df1

df1 %>%
  # mutate(facet = ifelse(name %in% c('HR11076', 'HR2476'), 'Low pH', 'Control')) %>%
  mutate(facet = ifelse(grepl("^HR110", name), '110 HPF', '24 HPF')) %>%
  mutate(facet = factor(facet, levels = c('24 HPF', '110 HPF'))) %>%
  mutate(moduleTraitCor = round(moduleTraitCor, 2)) %>%
  mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), moduleTraitCor)) %>%
  # mutate(star = ifelse(star != '', paste0(moduleTraitCor, '(', star,')'), '')) %>%
  ggplot(aes(y = module, x = name, fill = moduleTraitCor)) +
  geom_tile(color = 'white', size = 0.7, width = 1) +
  # geom_raster() +
  geom_text(aes(label = star),  vjust = 0.5, hjust = 0.5, size= 4, family =  "GillSans") +
  ggsci::scale_fill_gsea(name = "", reverse = T, na.value = "white") +
  # scale_fill_viridis_c(name = "Membership", na.value = "white") +
  ggh4x::scale_y_dendrogram(hclust = hclust) +
  labs(x = '', y = 'Module') +
  guides(fill = guide_colorbar(barwidth = unit(3.5, "in"),
    barheight = unit(0.1, "in"), label.position = "top",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.5,
    frame.colour = "black", frame.linewidth = 0.5,
    label.theme = element_text(size = 10))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white')) +
  facet_wrap(~ facet, scales = 'free_x') -> p1 

p1 <- p1 + theme(
  axis.line.x = element_blank(),
  axis.line.y = element_blank(),
  axis.text.y = element_text(hjust = 1.2),
  axis.ticks.length = unit(5, "pt"))

p1 <- p1 + theme(panel.spacing.x = unit(-0.5, "mm"))

# BARPLOT

p2 <- stats %>% 
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) %>%
  ggplot(aes(y = module)) +
  scale_x_continuous("Número de miRs") +
  geom_col(aes(x = n), width = 0.95, position = position_stack(reverse = TRUE)) +
  # geom_col(aes(x = reads_frac), width = 0.95, fill = "grey")
  # scale_fill_manual(name = '', values = c("#303960", "#647687", "#E7DFD5")) + # grey90
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    axis.title.y = element_blank(), 
    axis.text.y= element_blank(),
    axis.ticks.y =element_blank(), 
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.length = unit(5, "pt"))


library(patchwork)

# p1 + plot_spacer() + p2 + plot_layout(widths = c(5,-0.5, 10))  #& theme(plot.margin = 0)

psave <- p1 +  plot_spacer() + p2 + plot_layout(widths = c(7, -0.9, 3)) + labs(caption = '* corPvalueStudent < 0.05 ') 


ggsave(psave, filename = 'WGCNA_MIRS2HEATMAP.png', path = path, width = 7, height = 5, device = png, dpi = 300)

# network ====

module_members_1 <- c("grey", "black", "yellow", "pink", "turquoise")

module_members_2 <- c("brown", "red", "green", "blue")

library(igraph)
library(tidygraph)
library(ggraph)

exportNet <- function(TOMobj, moduleColors, threshold = 0.9) {
  
  nodeNames <- names(moduleColors)
  nodeAttr <- moduleColors
  
  
  # file_out <- gsub('.RData', '', basename(TOMFile))
  
  file_out <- "WGCNA.mirs"
    
  edgeFile <- paste0(path, "/",file_out, ".edges.txt")
  nodeFile <- paste0(path, "/",file_out, ".nodes.txt")
  
  
  Net = exportNetworkToCytoscape(TOMobj,
    edgeFile = edgeFile,
    nodeFile = nodeFile,
    weighted = TRUE,
    threshold = threshold,
    nodeNames = nodeNames, nodeAttr = nodeAttr)
  
  cat('\nEdges: ',nrow(Net$edgeData), 'and Nodes: ',nrow(Net$nodeData),'\n')
  
  graph <- tidygraph::tbl_graph(nodes = Net$nodeData, edges = Net$edgeData, directed = FALSE)
  
  return(graph)
  
  # if(is.null(query)) {
  #   graph
  # } else {
  #   graph %>% activate(nodes) %>% filter(nodeName %in% query)
  # }
  
  
  # return(Net)
  
}

g <- exportNet(TOM, moduleColors, threshold = 0.3)

str(SORT_MIRS <- c("MIR-278","LET-7","MIR-133",
  "Cluster_55760","Cluster_55776","MIR-2","MIR-315", "MIR-153", "BANTAM", "MIR-190", "MIR-2722", "MIR-1988", 
  "MIR-92", "MIR-277B", "MIR-216"))

# g %>% activate("nodes") %>% mutate(nodeName %in% SORT_MIRS)

scale_col <- g %>% activate("nodes") %>% as_tibble() %>% distinct(`nodeAttr[nodesPresent, ]`) %>% pull()

g <- g %>% activate("nodes") %>%  
  mutate(betweenness = betweenness(.), degree = centrality_degree(),
  membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
  pageRank = page_rank(.)$vector)

Levels <- c(module_members_2, module_members_1)

g <- g %>% activate("nodes") %>%
  mutate(`nodeAttr[nodesPresent, ]` = factor(`nodeAttr[nodesPresent, ]`, levels = Levels))

layout = create_layout(g, layout = 'igraph', algorithm = 'kk')


ggraph(layout) +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  geom_node_point(aes(color = `nodeAttr[nodesPresent, ]`, size = degree)) +
  # ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans") +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  theme(legend.position = "none") +
  coord_fixed() -> psleft

# ggsave(psave, filename = 'WGCNA_MIRS2NETWORK.png', path = path, width = 10, height = 10, device = png, dpi = 300)


node_labs <- paste(LETTERS[1:length(Levels)], ") ", Levels, sep = "")
names(node_labs) <- Levels


ggraph(layout) +
  facet_nodes(~ `nodeAttr[nodesPresent, ]`, labeller = labeller(`nodeAttr[nodesPresent, ]` = node_labs), 
    scales = "free") +
  scale_color_manual('', values = structure(scale_col, names = scale_col) ) +
  geom_edge_arc(aes(edge_alpha = weight), strength = 0.1) + # edge_width
  geom_node_point(aes(color = `nodeAttr[nodesPresent, ]`, size = degree)) +
  ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans") +
  scale_edge_width(range = c(0.3, 1)) +
  theme_graph(base_family = "GillSans") +
  guides(fill=guide_legend(nrow = 2)) +
  # coord_fixed() +
  guides(color = "none") +
  ggforce::geom_mark_hull(
    aes(x, y, group = as.factor(`nodeAttr[nodesPresent, ]`)),
    color = NA, fill = "grey76",
    concavity = 4,
    con.size = 0.3,
    con.linetype = 2,
    expand = unit(2, "mm"),
    alpha = 0.25)  +
  guides(fill = FALSE) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white')) -> ps


#
library(patchwork)

psave <- psleft + ps + patchwork::plot_layout(widths = c(0.7, 1.2))

ggsave(psave, filename = 'WGCNA_MIRS2NETWORK_2.png', path = path, width = 12, height = 8, device = png, dpi = 300)


# EXIT

g %>% 
  activate("nodes") %>% 
  as_data_frame(., "vertices") %>% 
  pivot_longer(cols = c('degree', 'betweenness')) %>%
  ggplot(aes(value, pageRank, color = as.factor(membership))) + 
  geom_point() +
  facet_grid(~name, scales = 'free')


g %>% 
  activate("nodes") %>% 
  as_data_frame(., "vertices") %>% 
  pivot_longer(cols = c('degree', 'betweenness', 'pageRank')) %>%
  ggplot(aes(value)) +
  facet_grid(~ name, scales = "free") +
  geom_histogram()

# separate go analysis by modules ====
# GO TO 5_SRNA_REGULATORY_FUNC.R

