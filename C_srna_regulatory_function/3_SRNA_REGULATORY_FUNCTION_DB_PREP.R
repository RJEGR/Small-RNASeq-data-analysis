
# BIND STEP 2_ W/ STEP 1)
# LOAD GENE ONTOLOGIES FROM BLAST (DIAMONT) BASED ON THE GENOMIC CDS.ALL.FA 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

# f <- "Trinotate_gene_ontology"

f <- "Trinotate_trans_ontology"

f <- list.files(path = paste0(wd, "ANNOTATIONS/"), pattern = f, full.names = T)

readMappings <-  function (file, sep = "\t", IDsep = ",")  {
  
  
  # From topGO::readMappings
  
  a <- read.delim(file = file, header = FALSE, quote = "", 
    sep = sep, colClasses = "character")
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1])
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, 
    split = IDsep)[[1]])))
}

MAP <- readMappings(f)

length(MAP) # 32672

genome_features <- read_rds(paste0(wd, "genome_features.rds"))

print(transcript_features <- genome_features[[1]])

print(gene_features <- genome_features[[2]])

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

nrow(out <- read_rds(paste0(wd, "SRNA_FUNCTION_PREDICTED.rds")))

# bind seqnames, gene_coords, gene_id, description, type, biotype to predicted srna regulatory func.

nrow(out <- out %>% left_join(gene_features)) # 6372

# bind gene ontology to transcript_id

# transcript_features <- transcript_features %>% mutate(transcript_id = sapply(strsplit(transcript_id, "[.]"), `[`, 1) )

str(transcript_features %>% distinct(transcript_id) %>% pull(.) -> query.ids)

n <- sum(keep <- names(MAP) %in% query.ids) / length(query.ids) 

cat('\n % of ids mapped in Gen Ontology: ', n*100, '%\n')

# % of ids mapped in Gen Ontology:  58.75308 %

paste_go <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}

TRANS2GO <- data.frame(transcript_id = rep(names(MAP[keep]),
  sapply(MAP[keep], length)),
  GO.ID = unlist(MAP[keep]), row.names = NULL) %>% 
  as_tibble() %>%
  group_by(transcript_id) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    n = n(),
    .groups = "drop_last")

transcript_features <- transcript_features %>% right_join(TRANS2GO)  
  
go_features <- transcript_features %>% dplyr::select(gene_id, GO.ID, n) %>% 
  group_by(gene_id) %>%
  summarise(
    n = sum(n),
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")

# bind go to predicted srna regulatory func.

sum(sort(unique(out$gene_id)) %in% sort(unique(go_features$gene_id))) # 3142 (58.75308 %)

out <- go_features %>% right_join(out, by = "gene_id")



# view(out)

# group go by query (i.e. by srna to prep for enrichment analysis)
# Here we look 147 miRs are well annotated (GaD)
# 

SRNA2GO <- out %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>% 
  distinct(query, GO.ID) %>%
  group_by(query) %>%
  summarise(
    across(GO.ID, .fns = paste_go), 
    .groups = "drop_last")


write_tsv(out, file = paste0(wd, "SRNA_REGULATORY_FUNCTION_DB.tsv"))

write_tsv(SRNA2GO, file = paste0(wd, "SRNA2GO.tsv"))

# GO TO DATAVIZ ====

nrow(SRNA2GO)

SRNA2GO <- split(strsplit(SRNA2GO$GO.ID, ";") , SRNA2GO$query)

SRNA2GO <- lapply(SRNA2GO, unlist)

typeof(SRNA2GO)

head(str(SRNA2GO))

# GO.ID ROWS IN SRNA 2 GO MUST BE A VECTOR AS IN MAP


typeof(MAP)

url <- "https://raw.githubusercontent.com/RJEGR/Cancer_sete_T_assembly/main/functions.R"

source(url)

# hsGO <- GOSemSim::godata('org.Hs.eg.db', ont="BP")

# write_rds(hsGO, paste0(wd, '/hsGO_BP.rds'))

hsGO <- read_rds(paste0(wd, '/hsGO_BP.rds'))

str(names(SRNA2GO) -> query.names)

query.p <- c(rep(0.05, length(query.names)))

allRes <- GOenrichment(query.p, query.names, SRNA2GO, Nodes = 20)

allRes %>% as_tibble() %>%
  arrange(Annotated) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  ggplot(aes(x = Term, y = Annotated, fill = p.adj.ks)) + # , fill = p.adj.ks
  coord_flip() +
  geom_col() +
  theme_minimal(base_family = "GillSans", base_size = 12)

# NOW SPLIT SRNAS BY DE AND INCLUDE TREATMENT TO TREATMENT GOENRICHMENT:
# EX.

SRNA2GO_ARM <- SRNA2GO %>% separate(query, into = c("query", "arm"), sep = "[.]")

SRNA2GO_ARM %>% distinct(arm) %>% pull() -> ARM

DF <- list()

for(i in ARM) {
  j <- i
  
  DF2GO <- SRNA2GO_ARM %>% filter(arm %in% j)
  
  nrow(DF2GO)
  
  DF2GO <- split(strsplit(DF2GO$GO.ID, ";") , DF2GO$query)
  
  DF2GO <- lapply(DF2GO, unlist)
  
  str(names(DF2GO) -> query.names)
  
  query.p <- c(rep(0.05, length(query.names)))
  
  allRes <- GOenrichment(query.p, query.names, DF2GO, Nodes = 20)
  
  
  allRes <- allRes %>% mutate(arm = j) %>% as_tibble()
  
  DF[[j]] <- allRes
  
}

do.call(rbind, DF) -> allRes

allRes %>% as_tibble() %>%
  arrange(Annotated) %>%
  mutate(Term = factor(Term, levels= unique(Term))) %>%
  ggplot(aes(x = Term, y = Annotated, fill = -log10(p.adj.ks))) + # , fill = p.adj.ks
  coord_flip() +
  geom_col() +
  theme_bw(base_family = "GillSans", base_size = 12) +
  facet_grid( arm ~ ., scales = "free") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ggsci::scale_fill_material()


allRes %>% 
  arrange(Annotated) %>%
  # mutate(Term = factor(Term, levels= unique(Term))) %>%
  # group_by(module, parentTerm) %>%
  # summarise(Annotated = sum(Annotated)) %>%
  # group_by(module) %>%
  # mutate(Size = Annotated / max(Annotated)) %>%
  mutate(Term = fct_reorder(Term, Annotated, .desc = F)) %>%
  ggplot(aes(y = Term, x = arm, size = Annotated)) +
  geom_point() + 
  theme_bw(base_family = "GillSans") +
  theme_bw(base_family = "GillSans", base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey86', color = 'white'),
    legend.position = 'top',
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ggsci::scale_fill_material()
