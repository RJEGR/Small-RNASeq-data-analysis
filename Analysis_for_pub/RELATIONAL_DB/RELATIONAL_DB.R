

# CREATE A RELATIONAL DB FOR SRNA_FUNCTIONAL AND RNA_LOCATION DATABASES == done
# VIEW IF INTRAGENIC ARE CO-REGULATING HOST GENE == zero
# VIEW IF TARGET IS COMMONLY EXPRESSED IN LARVAE DEV. == keep_expressed
# CORRELATE AND VISUALIZE MIR AND GENE EXPRESSION (MUST MATCH NEGATIVE CORRELATION)
# Graphical view of genomic location and regulation in genome = complicated!

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

pattern <- "multi_genome.newid.gff3$" # <- to get sequence-region per chromosome

wd_genome <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

# f <- list.files(path = wd_genome, pattern = pattern, full.names = T)

# G <- rtracklayer::import(f)

# length(.Genome <- G[G$type == "region"])

library(tidyverse)

# GR <- .Genome %>% as_tibble() %>% 
#   mutate(Chrom = gsub("^JALGQA010000", "", as.character(seqnames))) %>%
#   mutate(Chrom = gsub(".1$", "", Chrom)) %>%
#   dplyr::select(Chrom, start,  end, strand)
 
# GenomicRanges::GRanges(GR)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"


print(FUNCTIONAL_DB <- read_rds(paste0(wd,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))

is_na <- function(x) ifelse(is.na(x), 0, x)

keep_expressed <- FUNCTIONAL_DB %>% 
  dplyr::select(starts_with("SRR")) %>% 
  mutate(across(where(is.double), ~is_na(.))) %>% 
  rowSums()

FUNCTIONAL_DB <- FUNCTIONAL_DB[keep_expressed > 1,]


# Prepare miRNA_mRNA set -----

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

gene2tr <- read_rds(paste0(wd, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)


Lines <- FUNCTIONAL_DB %>% distinct(gene_id) %>%
  left_join(gene2tr) %>%   distinct(transcript_id) %>% pull()

writeLines(Lines, paste0(wd, "miRNA-mRNA.lines"))

wd <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# read_tsv(paste0(wd, "/RNA_LOCATION_DB.tsv")) %>% filter(SRNAtype == "miR")

# bind to Loci mirna db

print(.LOCATION_DB <- read_rds(paste0(wd, "RNA_LOCATION_MIR_DB.rds")))

MajorRNA2Name <- .LOCATION_DB %>% distinct(MajorRNA, Name) %>%   dplyr::rename("query" = "Name")

LOCATION_DB <- .LOCATION_DB %>%
  select(Name, MajorRNA, Locus, Reads, biotype_best_rank, gene_id) 

# 0 calculate distance between loci mirs ----

GenomicRanges::GRanges(LOCATION_DB)

# 1) VIEW IF INTRAGENIC ARE CO-REGULATING HOST GENE =====

sum(FUNCTIONAL_DB$gene_id %in% LOCATION_DB$gene_id)

FUNCTIONAL_DB %>%
  dplyr::select(!dplyr::starts_with("SRR")) %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  filter(arm == "mature") # must match only mature
  
.COORDS_TARGET <- FUNCTIONAL_DB %>%
  select(gene_id, gene_coords, query) %>%
  mutate(biotype_best_rank = "Target") %>%
  separate(gene_coords, into = c("Chrom", "Coords"), sep = "_") %>%
  separate(Coords, into = c("Start", "End", "Strand"), sep = ":") %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  separate(query, into = c("query", "arm"), sep = "[.]") %>%
  mutate(Chrom = gsub("^JALGQA010000", "", as.character(Chrom))) %>%
  mutate(Chrom = gsub(".1$", "", Chrom)) %>%
  select(-arm)

.COORDS_TARGET <- .COORDS_TARGET %>% left_join(MajorRNA2Name)

# a este punto no importa mucho el query

COORDS_TARGET <- .COORDS_TARGET %>%   select(-query) %>% distinct()
  
df <- .LOCATION_DB %>%
  dplyr::rename("query" = "Name") %>%
  select(any_of(names(COORDS_TARGET))) %>%
  rbind(COORDS_TARGET)

.LOCATION_DB %>% count(biotype_best_rank, sort = T)

df <- df %>% mutate_at(c("Start", "End"), as.numeric)



# from upgraded 
# find loci

df %>% 
  dplyr::filter(biotype_best_rank!= "Target") %>%
  distinct(MajorRNA, biotype_best_rank) %>%
  right_join(RES.P %>% distinct(MajorRNA, CONTRAST)) %>%
  group_by(CONTRAST) %>%
  count(biotype_best_rank)


# find target

# exit

data <- filter(GR, Chrom %in% unique(df$Chrom))

ggplot() +
  theme_classic(base_family = "GillSans", base_size = 12) +
  # geom_segment(data = df, aes(x = Start, xend = End, yend = Chrom, y = Chrom, color = biotype_best_rank), linewidth = 1.5, linetype="dashed") +
  geom_segment(data = data, aes(x = start, xend = end, yend = Chrom, y = Chrom), linewidth = 0.5, color="#E7DFD5", linetype="dashed") +
  gggenes::geom_gene_arrow(data = df, aes(xmin = Start, xmax = End, y = Chrom, fill = biotype_best_rank, color = biotype_best_rank)) 
  # facet_grid(~ biotype_best_rank)
  # gggenes::theme_genes()

# Create the chord plot (complicated)

library(GenomicRanges)

GR0 <- GRanges(Rle(df$Chrom), 
  ranges =  IRanges(start = df$Start, end = df$End), 
  strand = df$Strand, 
  biotype_best_rank = df$biotype_best_rank,
  gene_id = df$gene_id)

distances <- distanceToNearest(GR0)


df <- COORDS_TARGET %>% mutate_at(c("Start", "End"), as.numeric)

GR1 <- GRanges(Rle(df$Chrom), 
  ranges =  IRanges(start = df$Start, end = df$End), 
  strand = df$Strand, 
  biotype_best_rank = df$biotype_best_rank,
  gene_id = df$gene_id)

df <- .LOCATION_DB %>%
  dplyr::rename("query" = "Name") %>%
  select(any_of(names(COORDS_TARGET)))

GR2 <- GRanges(Rle(df$Chrom), 
  ranges =  IRanges(start = df$Start, end = df$End), 
  strand = df$Strand, 
  biotype_best_rank = df$biotype_best_rank,
  gene_id = df$gene_id)
# GenomicRanges::reduce(GR)

library(GenomicRanges)
z <- matrix(0,length(GR),3)
z[subjectHits(findOverlaps(GR1, GR)),1] <- 1
z[subjectHits(findOverlaps(GR2, GR)),2] <- 1
colnames(z) <- c("Genome", "Target", "Source")
mcols(GR) <- z


z <- matrix(0,length(GR2),3)
z[subjectHits(findOverlaps(GR1,GR2)),1] <- 1
z[subjectHits(findOverlaps(GR, GR2)),2] <- 1
colnames(z) <- c("mirnome", "Target", "Genome")

circlize::chordDiagramFromMatrix(z)  


library(ggforce)

ggplot(df) +
  geom_link(aes(x = Start, y = Chrom, xend = End, yend = Chrom, colour = biotype_best_rank, size = after_stat(indsex))) +
  ggforce::geom_arc(aes(x0 = 0, y0 = 0, start = Start, end = Chrom, r = 1))

# ggplot(df, aes(x = 0, y = 0, start_angle = Start, end_angle = End) 
  

# as network (olvida esto)

net <- .LOCATION_DB %>%
  dplyr::rename("query" = "Name") %>%
  select(any_of(names(.COORDS_TARGET))) %>%
  rbind(.COORDS_TARGET) %>%
  mutate(gene_id = ifelse(is.na(gene_id), biotype_best_rank, gene_id))

net <- net %>% distinct(gene_id, query, biotype_best_rank)

Node1 <- net %>% distinct(query, biotype_best_rank) %>% 
  dplyr::rename("Node" = "query")

Node2 <- net %>% distinct(gene_id, biotype_best_rank) %>% 
  dplyr::rename("Node" = "gene_id")


Nodes <- rbind(Node1,Node2)

Edges <- data.frame(from = net$gene_id, to = net$query)

graph <- tidygraph::tbl_graph(nodes = Nodes, edges = Edges, directed = T)

plot(graph)

library(igraph)
library(tidygraph)
library(ggraph)

# g <- graph %>% activate("nodes") %>%  
#   mutate(betweenness = betweenness(.), degree = centrality_degree(),
#     membership = igraph::cluster_louvain(., igraph::E(g)$weight)$membership,
#     pageRank = page_rank(.)$vector)

graph %>% activate("nodes") %>% filter(Node == "Cluster_18815")

layout <- create_layout(graph, layout = 'igraph', algorithm = 'kk')


ggraph(layout) +
  geom_edge_arc() +
  geom_node_point(aes(color = biotype_best_rank))
