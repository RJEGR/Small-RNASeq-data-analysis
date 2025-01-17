# TRACE Trinotate xls from transcriptome assembly
# Merge with Count-expression
# FIND transcrip/geneid in target analysis (RNAHybrid OR targetscan)
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

path_out <- wd

print(out <- read_tsv(paste0(wd, "SRNA_REGULATORY_FUNCTION_LONG_DB.tsv")))

# out <- out %>% filter(predicted == "BOTH")

# str(query.ids <- out %>% distinct(gene_id) %>% pull())


GTF <- "transcripts.gtf"

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/"

f <- list.files(path = wd, pattern = "METADATA_RNASEQ", full.names = T)

.colData <- read_csv(f)

f <- list.files(path = wd, pattern = GTF, full.names = T)

GTF <- rtracklayer::import(f)

print(GTF2DF <- GTF %>% 
    as_tibble() %>% 
    distinct(gene_id, transcript_id, ref_gene_id) %>% 
    filter(!is.na(ref_gene_id))) # A tibble: 70,001 × 3

# nrow(GTF2DF <- GTF2DF %>% filter(ref_gene_id %in% query.ids)) # 17532

# TRANSCRIPT ARE UNIQUE IDENTIFIER:
# gene_id ID AND ref_gene_id CLUSTERS TRANSCRIPT IDENTIFERS (SAME INFO)
# gene_id == ref_gene_id, EX. MSTRG.44204 == LOC124140464 W/ 41 DIFF. transcript_id (ISOFORMS)
# 

dir <- "~/Documents/MIRNA_HALIOTIS/"

wd <- paste0(dir, "RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/")

Manifest <- read_csv(list.files(path = wd, pattern = "METADATA_RNASEQ.csv", full.names = T))

Manifest[is.na(Manifest)] <- "Low CO2"

Manifest <- Manifest[!grepl("Mantle", Manifest$Isolated),]

Manifest <- mutate_if(Manifest, is.character, as.factor)



wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/QUANTIFICATION/"

f <- list.files(path = wd, pattern = "gene_count_matrix.csv", full.names = T)

.COUNT <- read_csv(f)

Mantle_sam <- c("SRR8956768", "SRR8956769")

# remove mantle samples

.COUNT <- .COUNT %>% dplyr::select(!starts_with(Mantle_sam)) 

is_na <- function(x) ifelse(is.na(x), 0, x)

keep_expressed <- .COUNT %>% 
  dplyr::select(starts_with("SRR")) %>% 
  mutate(across(where(is.double), ~is_na(.))) %>% 
  rowSums()

sum(keep_expressed <- keep_expressed > 1) # 131

nrow(.COUNT <- .COUNT[keep_expressed,]) 

.COUNT %>% distinct(gene_id) %>% nrow() # 58592 genes assembled. ie. expressed

any(GTF2DF$gene_id %in% .COUNT$gene_id)

GTF2DF <- GTF2DF[GTF2DF$gene_id %in% .COUNT$gene_id,]

COUNT <- GTF2DF %>% 
  distinct(gene_id, ref_gene_id) %>%
  right_join(.COUNT) %>%
  rename("assembled_id" = "gene_id", "gene_id" = "ref_gene_id")

which_cols <- COUNT %>% dplyr::select(dplyr::starts_with("SRR")) %>% names()

# optional, change transcript(ie. isoform ) to gene level

COUNT <- COUNT %>%
  group_by(gene_id) %>%
  summarise_at(vars(all_of(which_cols)), sum)

COUNT <- COUNT %>% drop_na(gene_id)

any(colnames(COUNT) %in% Manifest$LIBRARY_ID) # sanity check

# ANNOTATION ----

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/"

gene2tr <- read_rds(paste0(wd, "genome_features.rds"))[[1]] %>% 
  distinct(gene_id, transcript_id)

url <- "https://raw.githubusercontent.com/RJEGR/Small-RNASeq-data-analysis/master/FUNCTIONS.R"

source(url)

ref_path <- "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/ANNOTATIONS/trino"

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

go_f <- list.files(path = ref_path, pattern = "Trinotate_trans_ontology", full.names = T)

readMappings <-  function (file, sep = "\t", IDsep = ",")  {
  
  
  # From topGO::readMappings
  
  a <- read.delim(file = file, header = FALSE, quote = "", 
    sep = sep, colClasses = "character")
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1])
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, 
    split = IDsep)[[1]])))
}

bind_ids <- function(x) {
  x <- unique(x)
  x <- paste(x, sep = ';', collapse = ';')
}

MAP <- readMappings(go_f)

STRG2GO <- data.frame(transcript_id = rep(names(MAP),
  sapply(MAP, length)),
  GO.ID = unlist(MAP), row.names = NULL) %>% as_tibble() %>%
  left_join(gene2tr) %>%
  group_by(gene_id) %>%
  summarise(across(GO.ID, .fns = bind_ids)) 

.annot <- read_tsv(annot_f, na = ".") %>%
  distinct(transcript_id, sprot_Top_BLASTX_hit) # gene_ontology_BLASTX, Pfam, Kegg

.annot <- split_blast(.annot, "BLASTX")

annot <- .annot

# annot$`#gene_id` <- NULL


split_blast <- function (x, hit = "BLASTP", which_vars = "transcript_id") {
  
  # Upgraded function form trinotateR package:
  
  require(tidyverse)
  
  
  hit <- paste("sprot_Top_", hit, "_hit", sep = "")
  
  which_vars <- c(hit, which_vars) # "gene_id","transcript_id", "prot_id"
  
  y <- x %>% select_at(vars(all_of(which_vars))) %>% drop_na(any_of(hit))
  
  z <- y %>% pull(hit)
  
  z <- strsplit(z, "`")
  
  n <- sapply(z, length)
  
  z <- strsplit(unlist(z), "\\^")
  
  col2 <- sapply(z, "[", 2)
  
  col2[is.na(col2)] <- "."
  
  if (any(sapply(z, "[", 1) != col2))
    print("WARNING: check different values in columns 1 and 2")
  
  NAME <- gsub("^RecName: Full=", "", sapply(z, "[", 6))
  NAME <- gsub("SubName: Full=", "", NAME)
  NAME <- gsub(";$", "", NAME)
  NAME <- gsub(" \\{[^}]+}", "", NAME)
  
  # gene_id <- rep(y$gene_id, n)
  
  transcript_id <- rep(y$transcript_id, n)
  
  # protein <- rep(gsub(".*\\|", "", y$prot_id), n)
  
  uniprot <- sapply(z, "[", 1)
  
  align <- sapply(z, "[", 3)
  
  identity <- as.numeric(gsub("%ID", "", sapply(z, "[", 4)))
  
  evalue <- as.numeric(gsub("E:", "", sapply(z, "[", 5)))
  
  domain <- gsub("; .*", "", sapply(z, "[", 7))
  
  lineage <- sapply(z, "[", 7)
  
  genus <- gsub(".*; ", "", sapply(z, "[", 7))
  
  x1 <- data.frame(
    # gene_id,
    transcript_id, 
    # protein , 
    uniprot, align, identity, evalue, 
    name = NAME, lineage, domain, genus, stringsAsFactors = FALSE)
  
  message(nrow(x1), " ", hit, " annotations")
  
  as_tibble(x1)
}

any(gene2tr$transcript_id %in% annot$transcript_id)

annot <- gene2tr %>% 
  distinct(gene_id, transcript_id) %>%
  right_join(annot) %>%
  select(-transcript_id) %>% 
  distinct() %>%
  drop_na(uniprot)
  # rename("assembled_id" = "gene_id", "gene_id" = "ref_gene_id")

annot <- annot %>% 
  # filter(grepl("HUMAN", uniprot)) %>%
  mutate(uniprot = gsub("_[A-Z]*$", "", uniprot)) %>% distinct(gene_id, uniprot)

any(annot$gene_id %in% STRG2GO$gene_id)

# Generate DB with gene_id, uniprot and gene_ontology

annot <- annot %>% 
  group_by(gene_id) %>%
  summarise(across(uniprot, .fns = bind_ids), .groups = "drop_last") %>%
  left_join(STRG2GO)


DF <- annot %>% right_join(COUNT, by = "gene_id") 

write_tsv(DF, file.path(path_out, "Transcriptome_relational_db.tsv"))

# QUERY HOMEBOX ----

# The homeodomain is a highly conserved 60‐amino‐acid protein domain that is encoded by the homeobox and is found in organisms as diverse as mammals, insects, plants and yeast. Homeodomains function as DNA binding domains and are found in many transcription factors that control development and cell fate decisions

# bind_ids <- function(x) {
#   x <- unique(x)
#   x <- paste(x, sep = '^', collapse = '^')
# }

QUERY <- c("SIX1A^SIX1^SIX1B^SIX2^SIX6^SIX3^SIX4^HM33^HM32^SIX5^HM34^DUXA^HBX4^MIX2^DUX4C^DU4L4^DU4L2^DU4L3^DU4L5^DU4L6^DU4L7^DUX4^MIX1^DUX1^DUX5^MIXL1^SEBOX^DUX3^OTX2A^OTX2B^OTX2^OTX5^OTX1B^HME2A^HME1B^HME2B^HME1A^HME2^HMEN^HME1^HMIN^HM16^SMOX2^HME30^LHX9^LHX2^EVX1^HOX3^EVX2^VAB7^HM12^HMA2^VENT1^HXA2^HXB2^HXD3A^HXB2A^HXD3^HXA2B^PDX1^HXB3^CEH60^HXB7A^HXB7B^HXA7^HXB7^HXB6A^HXD5^HXC8^HXA6^HXB6B^HLOX2^ABDA^ARA^CAUP^LIM7^LHX5^AWH^LHX1^GSC^GSCA^GSCB^GSC2^ALX1^ALX3^VSX1^ALX4^SMOX5^NKX25^TLX2^TLX3^TLX1^BARX2^NKX32^BARX1^AL^ARX^RX3^ALX^RXA^RXB^RX1^ARXH^DRGX^PHX2B^EMX1^EMX2^HM02^VAX2^NOTO^VAX1B^VAX2B^VAX2A^VAX1^NOT2^HM09^VNT1B^LIM4^LHX6^LHX8^LHX3^LMX1A^LMX1B^ISX^PRRX2^RAX2^RX2^LHX4^HM14^VSX2^HHEX^HHEXH^HM7X^GBX1^HM90^GBX2^RHXF1^SHOX2^SHOX^CEH28^OTP^OTBP^PHX2A^CEH17^SMOX3^HXB6^DLLH^DLX1A^DLX4B^DLX6A^DLX4A^DLX1^DLL1^DLX6^BOX5^DLX2^DLX3^DLL4^DLX4^DLX2A^DLX2B^DLX5^DLL3^TGIF2^AKR^TGIF1^TF2LX^TF2LY^PKNX2^PKNX1^ME3L1^MEIS3^HOX71^MSX2^HMGX7^HOX7P^MSX1^VAX1A^BARH1^HESX1^ESX1^ANF1^NOBOX^KNOX3^KNOX4^HBX3^KNOS6^KNAP1^KNAT1^KNAP2^KNOX8^HM05^SAX1^HXD12^NKX23^EMS^NKX26^HXA4^HMD1^ZFHX4^ZFHX3^ZFHX2^HM07^LIM6^MEOX2^MEOX1^HXC4A^HXC4^HXB4^HXB5^HXD4^HXD4A^HXB4A^PITX^PITX2^PITX1^PITX3^UNC30^PV1^BSH^ZAX^HM30^KOZA^BARH2^HM01^NKX21^NKX24^HMH2^HM24^LOX10^HNK2^NKX22^NX22A^VND^NKX28^ZEB2^CUX1^CUX2^CUT^HM39^HM31^BAX1B^HMX3^HMX3A^HMX3B^HMX1^HMX^HMX2^NKX61^HM19^PNX^HXC11^HXCBA^HXA9A^PHP3^HXD11^HXABA^ABDB^HXDBA^HXD9^HXC9^HXB9^HXABB^HXA9^HXA11^HXC13^HXCDA^HXA13^ZEB1^ZAG1^HM10^RX^MEIS1^MEI3A^MEI3B^MEIS2^HTH^ME3L2^UNC62^HBX9^OTX^PROP1^PRRX1^SMOX4^HM13^HXA3^HXA3A^HXB3A^BAGP^HM08^NKX62^NKX63^NKX31^MSX3^UNPG^HXA2A^CDX1^CDX2^CDX4^CDX^HMD2^PAL1^SLOU^MSXB^UNC4^DVE1^HXA1^HXB1^HXD1^HXA1A^HXB1A^HXB1B^DBX1A^DBX1B^DBX1^DBX2^HLX^HOX1^GSX1^GSX2^BRN3^MNX1^MKX^EXD^HM20^HM40^HXA5^HXA5A^HXB5A^HXB5B^KNOX6^KNOX7^KNOX2^KNOS2^KNOSB^DMBX1^DMX1A^DMX1B^OTX5A^OTX5B^CRX^ROUGH^HM17^MSXA^HD12^OTX1^OTX1A^OTXH^HM37^PROS^PROS1^PROX1^PROX2^ONEC^HM38^HM21^HD2^HXD10^HMB4^HXA10^HXD9A^HXC10^HXAAB^HXC9A^HXDCA^HXA4A^SCR^HME60^LIN39^UNC42^HMBX1^TTX3")
# 
# .annot %>%  
#   # distinct(uniprot, name) %>%
#   filter(grepl("homeobox|Homeobox",name)) %>%
#   mutate(uniprot = gsub("_[A-Z]*$", "", uniprot)) %>%
#   distinct(uniprot) %>%
#   summarise(across(uniprot, .fns = bind_ids), .groups = "drop_last") %>%
#   pull()
# 

HOXDF <- gene2tr %>% 
  distinct(gene_id, transcript_id) %>%
  right_join(.annot) %>% 
  filter(grepl("homeobox|Homeobox",name)) %>%
  # ggplot(aes(identity)) + geom_histogram()
  mutate(uniprot = gsub("_[A-Z]*$", "", uniprot)) %>%
  distinct(gene_id, uniprot, name)

QUERY <- HOXDF %>% distinct(gene_id) %>% pull()

# QUERY <- unlist(strsplit(QUERY, "\\^"))

# HOMEBOX IS EXPRESSED  -----


any(DF$gene_id %in% HOXDF$gene_id)


DB1 <- DF %>% 
  # mutate(uniprot = strsplit(uniprot, ";")) %>%
  # unnest(uniprot) %>%
  filter(gene_id %in% QUERY)


# IS TARGETED??

any(out$gene_id %in% HOXDF$gene_id)

# MIR-10-5p is a member miR-10 family, known as a conserved Hox-cluster regulator.

MIR10 <- c("Cluster_42842", "Cluster_10732", "Cluster_10737")

Qtarget <- out %>%
  mutate(query = strsplit(query, ";")) %>%
  unnest(query) %>%
  mutate(query = gsub(".mature$","", query)) %>%
  filter(query %in% MIR10) %>%
  distinct(gene_id, query, predicted) %>%
  filter(gene_id %in% QUERY) %>%
  pull(gene_id)

HOXDF %>% filter(gene_id %in% Qtarget) %>% view()
  
DB1 %>% 
  filter(gene_id %in% Qtarget) %>%
  pivot_longer(cols = names(COUNT)[-1], values_to = "Reads", names_to = "LIBRARY_ID") %>%
  filter(Reads >0) %>%
  left_join(Manifest) %>% view()



out %>% distinct(gene_id, predicted) %>% count(gene_id, predicted, sort = T)

HOXDF %>% 
  left_join(COUNT, by = "gene_id") %>% 
  drop_na() %>% 
  # pivot_longer(cols = names(COUNT)[-1], values_to = "Reads", names_to = "LIBRARY_ID") %>%
  # filter(Reads >0) %>%
  # left_join(Manifest) %>%
  left_join(distinct(out, gene_id, predicted), by = "gene_id") %>% 
  # drop_na(predicted) %>% 
  # distinct()
  write_tsv(file.path(path_out, "Transcriptome_relational_hoxdb.tsv"))
