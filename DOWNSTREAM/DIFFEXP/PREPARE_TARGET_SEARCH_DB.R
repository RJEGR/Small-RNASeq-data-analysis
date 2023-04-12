
# IDENTIFICAR PROTEINAS Y DOMINIOS PROTEICOS DE LA BIOMINERALIZACION
# 

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)


library(tidyverse)

wd <- "~/Documents/MIRNA_HALIOTIS/"

ref_path <- paste0(wd, "/RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/ANNOTATION/")

annot_f <- list.files(path = ref_path, pattern = "Trinotate.xls", full.names = T)

go_f <- list.files(path = ref_path, pattern = "Trinotate_report.xls.gene_ontology", full.names = T)

readMappings <-  function (file, sep = "\t", IDsep = ",")  {
  
  
  # From topGO::readMappings
  
  a <- read.delim(file = file, header = FALSE, quote = "", 
    sep = sep, colClasses = "character")
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1])
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, 
    split = IDsep)[[1]])))
}

# MAP <- readMappings(go_f)

annot <- read_tsv(annot_f, na = ".")

annot <- annot %>% drop_na(sprot_Top_BLASTP_hit)

names(annot)[1] <- "gene_id"

annot %>% head() %>% view()

blastp_df <- split_blast(annot, hit = "BLASTP")

go_df <- split_gene_ontology(annot, hit = "BLASTP")

# EXPLORATORY   ======

go_df %>% count(ontology)

hist(blastp_df$identity)

prot_list <- c("Perlucin", "Protein PIF", "Carbonic anhydrase", "Chitin synthase", "Insoluble matrix shell protein")

# prot_list <- str_to_lower(prot_list)

prot_list <- paste(prot_list, collapse = "|")

blastp_df %>% 
  # mutate_at(vars(all_of("name")), list(~ str_to_lower(.))) %>%
  filter_all(any_vars(grepl(prot_list, .))) %>% view() 
# mutate_at(vars(all_of(into)), list(~ str_to_sentence(.)))
# blastp_df %>% filter_all(any_vars(grepl("LOC", .))) %>%
  # distinct(gene, transcript, protein)

str(query_protein <- blastp_df %>% 
  filter_all(any_vars(grepl(prot_list, .))) %>%
  distinct(protein) %>% pull())


f <- list.files(path = ref_path, pattern = "transcripts.fa.transdecoder.rmdup.cds", full.names = T)


seqs <- Biostrings::readDNAStringSet(f)

str(keep <- names(seqs))

keep <- sapply(strsplit(keep, " "), `[`, 1)

str(keep <- unique(keep))

sum(keep <- keep %in% query_protein)

subseqs <- seqs[keep]# <--------

# blastp_df %>% distinct(name) %>% view()

blastp_df %>% distinct(uniprot) %>% view()

blastp_df %>% distinct(genus) %>% view()

blastp_df %>% 
  filter_all(any_vars(grepl("Mollusca", .))) %>%
  count(genus) %>% view()


blastp_df %>% 
  filter_all(any_vars(grepl("Mollusca", .))) %>% view()





# GROUP GENES BY THEMATIC COMPONENTS ======
# EX. KEGG, BIOLOGICAL PROCESSS, CELLULAR COMPONENTS OR MOLEC. FUNCTION.

# GENE SAMPLE ABUNDANCE  ======



