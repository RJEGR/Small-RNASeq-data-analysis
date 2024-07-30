# WGCNA FOR ONLY MIRS
# this version use the majorRNA sequences summarization
# 
# wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"
# read_tsv(paste0(wd, "DESEQ2SRNATARGET.tsv"))

library(WGCNA)
library(flashClust)
library(tidyverse)


rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

f <- file.path(path, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds")

datExpr <- read_rds(f)


# Nota: usar transformacion ablanda los datos
# datExpr <- log2(datExpr+1)
# z_scores <- function(x) {(x-mean(x))/sd(x)}
# datExpr <- apply(datExpr, 1, z_scores)

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

# power_pct <- quantile(soft_values, probs = 0.95)

power_pct <- soft_values[which.max(soft_values)]

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
  facet_grid(~name, scales = 'free_x', switch = "x") +
  geom_text(aes(label = Power), size = 2) +
  geom_abline(slope = 0, intercept = softPower, linetype="dashed", alpha=0.5) +
  # geom_vline(xintercept = min(meanK), linetype="dashed", alpha=0.5) +
  labs(y = 'Soft Threshold (power)', x = '', 
    caption = caption) +
  # scale_x_continuous(position = "top") +
  theme_light(base_family = "GillSans",base_size = 16) -> psave

psave

# Build a adjacency "correlation" matrix

enableWGCNAThreads()

# specify network type

adjacency <- adjacency(datExpr, 
  power = softPower, 
  type = "unsigned")


TOM <- TOMsimilarity(adjacency, TOMType = "unsigned") # specify network type


dissTOM = 1 - TOM

heatmap(TOM)

# Generate Modules ----
# Generate a clustered gene tree

geneTree = flashClust(as.dist(dissTOM), method="average")

plot(geneTree, xlab="", sub="", 
  main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#This sets the minimum number of genes to cluster into a module

#   Modules of co-expressed miRNAs were then determined by using the dynamic tree cut procedure with a minimum module size cut-off of 5. This cut-off was chosen considering the miRNA transcriptome’s small size and the fact that a single miRNA can target multiple RNA transcripts. With this procedure, 18 different modules of EigenmiRNA were identified (https://doi.org/10.1038/s41598-020-)

minClusterSize <- softPower

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

MEDissThres = 0

title <- paste0("Dynamic Tree Cut", "Merged dynamic\n(cutHeight=", MEDissThres,")")
  
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors = merge$colors

mergedMEs = merge$newMEs

#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, dynamicColors, c("Modules"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), title, dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors <- mergedColors

colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1


MEs <- mergedMEs


# how modules where obtained:
nm <- table(moduleColors)


cat('Number of mudules obtained\n :', length(nm))

print(nm)
MEDissThres

WGCNA <- as_tibble(moduleColors, rownames = "MajorRNA") %>% 
  dplyr::rename("WGCNA" = "value")

WGCNA <- read_tsv(list.files(path = path, 
  pattern = "SEQUENCES_MERGED_DESEQ_RES.tsv", full.names = T)) %>%
  left_join(WGCNA)

write_rds(WGCNA, file = paste0(path, "SEQUENCES_MERGED_DESEQ_RES_WGCNA.rds"))

# exit

# Plot TOM
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
TOMplot(plotTOM, dendro = geneTree, Colors = moduleColors)

#

# LABEL MIRS?

RES.P <- read_tsv(paste0(path, "SEQUENCES_MERGED_DESEQ_RES.tsv")) %>% 
  filter( padj < 0.05 & abs(log2FoldChange) > 1)

query.names_1 <- RES.P %>% filter(CONTRAST %in% "CONTRAST_A") %>% 
  filter(log2FoldChange < 0) %>%
  distinct(MajorRNA) %>% pull()

query.names_2_up <- RES.P %>% filter(CONTRAST %in% "CONTRAST_B") %>% 
  filter(log2FoldChange < 0) %>%
  distinct(MajorRNA) %>% pull()

query.names_2_down <- RES.P %>% filter(CONTRAST %in% "CONTRAST_B") %>% 
  filter(log2FoldChange > 0) %>%
  distinct(MajorRNA) %>% pull()

reads <- colSums(datExpr)

Total <- sum(reads)

# Sanity check:

identical(names(colSums(datExpr)), names(moduleColors))

data.frame(reads, moduleColors) %>% 
  as_tibble(rownames = "Name") %>% 
  mutate(DE = "Basal") %>%
  mutate(DE = ifelse(Name %in% query.names_1, "24 HPF (Up)", DE)) %>%
  mutate(DE = ifelse(Name %in% query.names_2_up, "110 HPF (Up)", DE)) %>%
  mutate(DE = ifelse(Name %in% query.names_2_down, "110 HPF (Down)", DE)) %>%
  group_by(moduleColors, DE) %>% 
  summarise(n = n(), reads = sum(reads)) %>%
  dplyr::rename('module' = 'moduleColors') -> stats

# INSTEAD OF BINARY USE DATA FROM RESPIROMETRY AND MORPHOLOGY

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"
# 
datTraits <- read_tsv(paste0(path, "/WGCNA_datTraits.tsv"), )

which_cols <- colnames(datTraits)[!grepl(c("^HR"), colnames(datTraits))]  

datTraits <- datTraits %>% dplyr::select_at(which_cols)

row_names <- datTraits$Row

identical(rownames(datExpr),row_names)

keep <- grepl(c("HR11076|HR2476"), row_names)


datTraits <- datTraits %>% filter(Row %in% row_names[keep])

datTraits$Row <- NULL

datTraits <- as(datTraits, "matrix")

rownames(datTraits) <- row_names[keep]

# rownames(datTraits) <- row_names


# OR Prep binary datTraits

# datTraits <- gsub(".clean.newid.subset", "", rownames(datExpr))
# 
# HR11076 <- grepl('HR11076', datTraits)
# 
# HR1108 <- grepl('HR1108', datTraits)
# 
# HR2476 <- grepl('HR2476', datTraits)
# 
# HR248 <- grepl('HR248', datTraits)
# 
# datTraits <- data.frame(HR11076, HR1108, HR2476, HR248)
# 
# datTraits <- 1*datTraits

# rownames(datTraits) <- rownames(datExpr)

# Recalculate MEs with color labels

# MEs0 = moduleEigengenes(datExpr[,keep], moduleColors)$eigengenes

MEs0 = moduleEigengenes(datExpr[keep,], moduleColors)$eigengenes

MEs = orderMEs(MEs0)

names(MEs) <- str_replace_all(names(MEs), '^ME', '')

# heatmap(as(MEs, "matrix"))

moduleTraitCor = WGCNA::cor(MEs, datTraits, use= "p")

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
    strip.background = element_rect(fill = 'grey89', color = 'white')) 
# facet_grid(~ facet, scales = 'free_x',space = "free_x")

# recode_to <- structure(c("Development", "Growth", "Calcification", "Respiration", "pH 7.6", "pH 8.0", "pH 7.6", "pH 8.0"), names = colnames(datTraits))

recode_to <- c("Desarrollo", "Crecimiento", "Calcificación", "Respiración")
recode_to <- structure(recode_to,names = colnames(datTraits))


df1 %>%
  mutate(facet = "") %>%
  mutate(facet = ifelse(!grepl("^HR", name), 'A) Evaluación', facet)) %>%
  mutate(facet = ifelse(grepl("^HR24", name), 'B) 24 HPF', facet)) %>% 
  mutate(facet = ifelse(grepl("^HR110", name), 'C) 110 HPF', facet)) %>% 
  mutate(facet = factor(facet, levels = c('A) Evaluación', 'B) 24 HPF', 'C) 110 HPF'))) %>%
  mutate(name = recode_factor(name, !!!recode_to, .ordered = T)) %>%
  # mutate(name = factor(name, levels = )) %>%
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
    strip.background = element_rect(fill = 'grey89', color = 'white')) +
  facet_grid(~ facet, scales = 'free_x',space = "free_x") -> p1

p1 <- p1 + theme(
  axis.line.x = element_blank(),
  axis.line.y = element_blank(),
  axis.text.y = element_text(hjust = 1.2),
  axis.ticks.length = unit(5, "pt"))

p1 <- p1 + theme(panel.spacing.x = unit(-0.5, "mm"))

# BARPLOT


p2 <- stats %>% 
  mutate(module = factor(module, levels = hclust$labels[hclust$order])) %>%
  ggplot(aes(y = module)) + #  fill = DE, color = DE
  scale_x_continuous("Número de miRs") +
  geom_col(aes(x = n), width = 0.95, position = position_stack(reverse = TRUE)) +
  # geom_col(aes(x = reads_frac), width = 0.95, fill = "grey")
  # scale_fill_manual(name = '', values = c("#303960", "#647687", "#E7DFD5")) + # grey90
  theme_classic(base_size = 14, base_family = "GillSans") +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.title.y = element_blank(), 
    axis.text.y= element_blank(),
    axis.ticks.y =element_blank(), 
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.length = unit(5, "pt"))


library(patchwork)

# p1 + plot_spacer() + p2 + plot_layout(widths = c(5,-0.5, 10))  #& theme(plot.margin = 0)

psave <- p1 +  plot_spacer() + p2 + plot_layout(widths = c(7, -0.25, 1.5)) + labs(caption = '* corPvalueStudent < 0.05 ') 

# psave 

ggsave(psave, filename = 'WGCNA_MIRS2HEATMAP_2.png', path = path, width = 10, height = 5, device = png, dpi = 300)

# network ====

# module_members_1 <- c("grey", "black", "yellow", "pink", "turquoise")

# module_members_2 <- c("brown", "red", "green", "blue")

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

str(SORT_MIRS <- c("MIR-278","LET-7","MIR-133",
  "Cluster_55760","Cluster_55776","MIR-2","MIR-315", "MIR-153", "BANTAM", "MIR-190", "MIR-2722", "MIR-1988", 
  "MIR-92", "MIR-277B", "MIR-216"))

# g <- exportNet(TOM, moduleColors, threshold = 0.3)

g <- exportNet(TOM, moduleColors, threshold = 0)

g <- g %>% activate("edges") %>% mutate(weight = ifelse(weight < 0.3, NA, weight)) %>% filter(!is.na(weight))

moduleColors[names(moduleColors) %in% SORT_MIRS] %>% as_tibble(rownames = "Family") %>% view()

# g1 <- g %>% activate("nodes") %>% 

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
  ggrepel::geom_text_repel(data = layout, aes(x, y, label = nodeName), max.overlaps = 50, family = "GillSans", size = 2) +
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

