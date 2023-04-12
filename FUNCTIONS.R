
# Evaluate: https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/functions.R

split_blast <- function (x, hit = "BLASTP") {
  
  # Upgrade function form trinotateR package:
  
  require(tidyverse)
  
  
  hit <- paste("sprot_Top_", hit, "_hit", sep = "")
  
  which_vars <- c(hit, "gene_id", "transcript_id", "prot_id")
  
  y <- x %>% select_at(vars(all_of(which_vars))) %>% drop_na(any_of(hit))
  
  z <- y %>% pull(hit)
  
  z <- strsplit(z, "`")

  n <- sapply(z, length)
  
  z <- strsplit(unlist(z), "\\^")
  
  if (any(sapply(z, "[", 1) != sapply(z, "[", 2)))
    print("WARNING: check different values in columns 1 and 2")
  
  NAME <- gsub("^RecName: Full=", "", sapply(z, "[", 6))
  NAME <- gsub("SubName: Full=", "", NAME)
  NAME <- gsub(";$", "", NAME)
  NAME <- gsub(" \\{[^}]+}", "", NAME)
  
  gene <- rep(y$gene_id, n)
  
  transcript <- rep(y$transcript_id, n)
  
  protein <- rep(gsub(".*\\|", "", y$prot_id), n)
  
  uniprot <- sapply(z, "[", 1)
  
  align <- sapply(z, "[", 3)
  
  identity <- as.numeric(gsub("%ID", "", sapply(z, "[", 4)))
  
  evalue <- as.numeric(gsub("E:", "", sapply(z, "[", 5)))
  
  domain <- gsub("; .*", "", sapply(z, "[", 7))
  
  lineage <- sapply(z, "[", 7)
  
  genus <- gsub(".*; ", "", sapply(z, "[", 7))
  
  x1 <- data.frame(gene, transcript, protein , 
                  uniprot, align, identity, evalue, 
                  name = NAME, lineage, domain, genus, stringsAsFactors = FALSE)
  
  message(nrow(x1), " ", hit, " annotations")
  
  as_tibble(x1)
}

split_gene_ontology <- function(x, hit = "BLASTP") {
  
  require(tidyverse)
  
  gene_ontology_hit <-  paste("gene_ontology_", hit, sep = "")
  
  which_vars <- c(gene_ontology_hit, "gene_id", "transcript_id", "prot_id")
  
  y <- x %>% select_at(vars(all_of(which_vars))) %>% drop_na(any_of(gene_ontology_hit))
  
  z <- y %>% pull(gene_ontology_hit)
  
  z <- strsplit(z, "`")
  
  n <- sapply(z, length)
  
  z <- strsplit(unlist(z), "\\^")
  
  x1 <- data.frame(gene = rep(y$gene_id, n), 
    transcript = rep(y$transcript_id, n), 
    protein = rep(gsub(".*\\|", "", y$prot_id), n), 
    go = sapply(z, "[", 1), ontology = sapply(z, "[", 2), 
    name = sapply(z, "[", 3), stringsAsFactors = FALSE)
  
  message(nrow(x1), " ", gene_ontology_hit, " annotations")
  
  as_tibble(x1)
  
}

split_KEGG <- function (x, hit = "Kegg")  {
  library(data.table)
  y <- x[!is.na(get(hit)), .(get(hit), gene_id, transcript_id, 
    prot_id)]
  z <- strsplit(y$V1, "`")
  n <- sapply(z, length)
  z <- strsplit(unlist(z), "\\^")
  x1 <- data.frame(gene = rep(y$gene_id, n), 
    transcript = rep(y$transcript_id, n), 
    protein = rep(gsub(".*\\|", "", y$prot_id), n), 
    Kegg = sapply(z, "[", 1), stringsAsFactors = FALSE)
  message(nrow(x1), " ", hit, " annotations")
  data.table(x1)
}
