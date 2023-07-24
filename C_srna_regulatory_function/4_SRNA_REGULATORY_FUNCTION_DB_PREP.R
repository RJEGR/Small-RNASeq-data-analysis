# JOIN TRANSCRIPTOMIC DATA RELATED TO MIR-TARGETS:
# /Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

GTF <- "transcripts.gtf"

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"

f <- list.files(path = wd, pattern = GTF, full.names = T)

GTF <- rtracklayer::import(f)

# GTF %>% as_tibble() %>% dplyr::count(type)

GTF <- GTF[GTF$type == "transcript"]

seqlevels(GTF) <- unique(unfactor(seqnames(GTF)@values))

# Remember:

length(GTF) # 149366

length(IRanges::ranges(GTF)) # 149366

length(IRanges::reduce(GTF)) # 56566

# which is close to :

length(sort(unique(GTF$gene_id )))

length(range(GTF)) # 568

# ES NECESARIO REDUCIR EL OBJ, GTF A UN SOLO IDENTIFICADOR (GEN) POR COORDENADA:

# GTF[GTF$gene_id == "MSTRG.15"]
IRanges::reduce(GTF[GTF$gene_id == "MSTRG.15"])

# IRanges::ranges(GTF[GTF$gene_id == "MSTRG.15"], use.mcols = T)

GTFr <- IRanges::reduce(GTF)

# IRanges::subsetByOverlaps(GTFr,  GTF)

# IRanges::mergeByOverlaps(GTF, GTFr)

# IRanges::resize(GTF, width(range(GTF)), fix = "start")[GTF$gene_id == "MSTRG.15"]

# BIND W/ ANNOTATION FIRST

load(paste0(wd, "ANNOTATION/OUTPUTS/annot.Rdata"))

names(blastp_df)[1:2] <- c("gene_id", "transcript_id")

blastp_spiralia <- blastp_df %>% filter(grepl("Spiralia", lineage)) 

blastp_df <- blastp_df %>% 
  anti_join(blastp_spiralia) %>%
  group_by(transcript_id) %>% 
  arrange(desc(identity)) %>% 
  sample_n(1) %>% arrange(desc(identity)) %>%
  rbind(blastp_spiralia) %>%
  select(gene_id, transcript_id, protein, uniprot, genus, name) 

any(sort(unique(blastp_df$gene_id)) %in% sort(unique(GTF$gene_id)))

# BECAUSE ISOFORM DUPLICATES:
# EX.
# GTF %>% as_tibble() %>%  count(gene_id, sort = T)
# GTF %>% as_tibble() %>% filter(gene_id == "MSTRG.15")

gene_features_orf_rnaseq <- GTFr %>% as_tibble() %>% 
  left_join(GTF %>% as_tibble() %>% 
      select(start, end, gene_id), by = c("start", "end")) %>%
  mutate(gene_coords = paste(start, end, strand, sep = ":")) %>%
  select(seqnames, gene_coords, gene_id) %>%
  filter(gene_id %in% sort(unique(blastp_df$gene_id))) %>%
  distinct(seqnames, gene_coords, gene_id) %>%
  left_join(blastp_df, by = "gene_id") %>% 
  distinct() 

write_rds(gene_features_orf_rnaseq, paste0(wd, "gene_features.rds"))

# THOSE COORDS MAP W/ MIR TARGET?

wd1 <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT"

df1 <- read_rds(paste0(wd1, "/mir_vs_utr_rmdup_RNAhybrid.out.psig.gene_features.rds"))

df1 %>% left_join(GTF %>% as_tibble())

subsetByOverlaps(x = GTF, ranges = assembly[6,], type="any") %>% 
  as_tibble() %>% 
  mutate(Name = "Cluster_6") %>%
  select(seqnames, start, end, strand, type, biotype, Name, gene_id, transcript_id) 


# INCLUDE GO.

f <- list.files(path = paste0(wd, "ANNOTATION/OUTPUTS/"), 
  pattern = "Trinotate_report.xls.gene_ontology", full.names = T)


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

