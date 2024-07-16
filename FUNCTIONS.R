
# Evaluate: https://github.com/RJEGR/Cancer_sete_T_assembly/blob/main/functions.R

split_blast <- function (x, hit = "BLASTP") {
  
  # Upgraded function form trinotateR package:
  
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
  
  # Upgraded function form trinotateR package:
  
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

split_pfam <- function(x, hit = "Pfam"){
  
  # Upgraded function form trinotateR package:
  
  require(tidyverse)
  
  which_vars <- c(hit, "gene_id", "transcript_id", "prot_id")
  
  y <- x %>% select_at(vars(all_of(which_vars))) %>% drop_na(any_of(hit))
  
  z <- y %>% pull(hit)
  z <- strsplit(z, "`")

  n <- sapply(z, length)
  
  # split annotation into 5 columns
  
  z <- strsplit( unlist(z), "\\^" )
  
  gene <- rep(y$gene_id, n)
  
  transcript <- rep(y$transcript_id, n)
  
  protein <- rep(gsub(".*\\|", "", y$prot_id), n)
  
  pfam <- gsub("\\.[0-9]*", "", sapply(z, "[", 1))
  
  x1 <- data.frame(gene, 
    transcript, protein, pfam, 
    symbol = sapply(z, "[", 2), 
    name = sapply(z, "[", 3), 
    align = sapply(z, "[", 4), 
    evalue = as.numeric(gsub("E:", "", sapply(z, "[", 5) )), 
    stringsAsFactors = FALSE)
  
  message(nrow(x1), " ", hit, " annotations")
  
  as_tibble(x1)
  
}

split_kegg <- function (x, hit = "Kegg")  {
  # library(data.table)
  
  # y <- x[!is.na(get(hit)), .(get(hit), gene_id, transcript_id, prot_id)]
  
  which_vars <- c(hit, "gene_id", "transcript_id", "prot_id")
  
  y <- x %>% select_at(vars(all_of(which_vars))) %>% drop_na(any_of(hit))
  
  z <- y %>% pull(hit)
  
  z <- strsplit(z, "`")
  
  n <- sapply(z, length)
  z <- strsplit(unlist(z), "\\^")
  x1 <- data.frame(gene = rep(y$gene_id, n), 
    transcript = rep(y$transcript_id, n), 
    protein = rep(gsub(".*\\|", "", y$prot_id), n), 
    Kegg = sapply(z, "[", 1), stringsAsFactors = FALSE)
  message(nrow(x1), " ", hit, " annotations")
  
  as_tibble(x1)
}

get_eggnog <- function (x, ids, by = "transcript_id")  {
  # trinotateR::plot_NOGs()
  
  if (by == "transcript_id") {
    
    x1 <- x %>% filter(transcript_id %in% ids)
    
    # y <- unique(x1[!is.na(eggnog), .(transcript_id, eggnog)])
    eggnog <- x1 %>% drop_na(eggnog) %>% distinct(eggnog) %>% pull(eggnog)
    
  }
  else {
    x1 <- x %>% filter(gene_id %in% ids)
    
    # y <- unique(x1[!is.na(eggnog), .(gene_id, eggnog)])
    
    eggnog <- x1 %>% drop_na(eggnog) %>% distinct(eggnog) %>% pull(eggnog)
    
  }
  
  nogs <- gsub("(.*)\\^.*", "\\1", eggnog)
  
  # nogs <- eggnog
  
  
  # y %>% separate(col = eggnog, sep = "(.*)\\^.*", into = c('nogs', 'eggnog'))
  
  n <- match(nogs, egg$nog)
  
  y <- table(unlist(strsplit(egg$class[n], "")))
  
  y <- data.frame(y)
  
  names(y) <- c('code', 'Freq')
  
  
  
  return(y)
}


runtopGO <- function(topGOdata, topNodes = 20, conservative = TRUE) {
  
  RFisher <- runTest(topGOdata, 
    algorithm = "classic", 
    statistic = "fisher")
  
  # To make this test conservative. Next we will test the enrichment using the Kolmogorov-Smirnov test. We will use the both the classic and the elim method.
  
  if(conservative) 
  {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    RKS.elim <- runTest(topGOdata, algorithm = "elim", 
      statistic = "ks")
    
    
    
    
    allRes <- GenTable(topGOdata, 
      classicFisher = RFisher,
      classicKS = RKS, 
      elimKS = RKS.elim,
      orderBy = "elimKS", 
      ranksOf = "classicFisher", 
      topNodes = topNodes) 
  } else {
    RKS <- runTest(topGOdata, algorithm = "classic", 
      statistic = "ks")
    
    test.stat <- new("weightCount",
      testStatistic = GOFisherTest,
      name = "Fisher test", sigRatio = "ratio")
    
    weights <- getSigGroups(topGOdata, test.stat)
    
    allRes <- GenTable(topGOdata,
      classic = RFisher,
      KS = RKS,
      weight = weights,
      orderBy = "weight",
      ranksOf = "classic",
      topNodes = topNodes)
    
    # allRes <- GenTable(object = topGOdata, 
    #                    elimFisher = RFisher,
    #                    topNodes = topNodes)
  }
  
  return(allRes)
}

GOenrichment <- function(query.p, query.names, gene2GO, cons = T, onto = "BP", Nodes = Inf) {
  
  require(topGO)
  
  names(query.p) <- query.names
  
  
  # keep MAP of query genes
  
  keep <- names(gene2GO) %in% names(query.p) 
  
  gene2GO <- gene2GO[keep]
  
  keep <- names(query.p) %in% names(gene2GO)
  
  query.p <- query.p[keep]
  
  description <- "complete topGO enrichment using split_annot"
  
  
  topGOdata <- new("topGOdata", 
    ontology = onto, 
    description = description,
    allGenes = query.p,
    # geneSel = function(x) { x == 1 },
    geneSel = function(x) x,
    annot = annFUN.gene2GO,
    # mapping = hsGO, # omit this flag
    gene2GO = gene2GO)
  
  # run TopGO results 
  
  allGO <- usedGO(topGOdata)
  
  if(is.infinite(Nodes)) {
    topNodes <- length(allGO)
  } else {
    topNodes <- Nodes
  }
  
  
  allRes <- runtopGO(topGOdata, topNodes = Nodes, conservative = cons)
  
  # make p adjustable
  
  p.adj.ks <- p.adjust(allRes$classicKS , method="BH")
  
  allRes <- cbind(allRes, p.adj.ks)
  
  allRes$Term <- gsub(" [a-z]*\\.\\.\\.$", "", allRes$Term)
  allRes$Term <- gsub("\\.\\.\\.$", "", allRes$Term)
  
  return(allRes)
  
  
}

get_res <- function(dds, contrast, alpha_cutoff = 0.1) {
  
  sA <- contrast[1]
  sB <- contrast[2]
  
  contrast <- as.character(DESeq2::design(dds))[2]
  
  keepA <- as.data.frame(colData(dds))[,contrast] == sA
  keepB <- as.data.frame(colData(dds))[,contrast] == sB
  
  contrast <- c(contrast, sA, sB)
  
  res = results(dds, contrast, alpha = alpha_cutoff)
  
  
  baseMeanA <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepA])
  baseMeanB <- rowMeans(DESeq2::counts(dds,normalized=TRUE)[,keepB])
  
  res %>%
    as.data.frame(.) %>%
    cbind(baseMeanA, baseMeanB, .) %>%
    cbind(sampleA = sA, sampleB = sB, .) %>%
    as_tibble(rownames = "Name") %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    mutate_at(vars(!matches("Name|sample|pvalue|padj")),
      round ,digits = 2)
}

prep_DE_data <- function(res, alpha, lfcThreshold) {
  
  FC_cols <- c('logFC','log2FoldChange')
  pv_cols <- c('P.Value','pvalue', 'PValue')
  pvad_cols <- c('padj', 'FDR', 'adj.P.Val')
  
  sam <- c('sampleA',  'sampleB')
  samv <- c('baseMeanA', 'baseMeanB')
  
  rename_to <- c('logFC', 'pvalue', 'padj')
  # 
  # sigfc <- "<b>P</b>-value & Log<sub>2</sub> FC"
  # pv <- "<b>P</b>-value"
  # fc <- "Log<sub>2</sub> FC"
  # 
  # sigfc <- expression(p - value ~ and ~ log[2] ~ FC)
  # pv <- "p-value"
  # fc <- expression(Log[2] ~ FC)
  # 
  sigfc <- "p - value ~ and ~ log[2] ~ FC"
  pv <- "p-value"
  fc <- "Log[2] ~ FC"
  
  
  sample_NS <- function(res, n) {
    
    x <- filter(res, !(cc %in% c('NS', fc)))
    
    y <- sample_n(filter(res, cc == fc), n)
    
    z <- sample_n(filter(res, cc == 'NS'), n) 
    
    dat_sampled <- rbind(x,y,z)
    
    return(dat_sampled)
  }
  
  res %>%
    # drop_na() %>%
    # select_at(vars(Name, contains(sam), contains(samv), contains(FC_cols), contains(pv_cols), contains(pvad_cols))) %>%
    rename_at(vars(contains(FC_cols),contains(pv_cols),contains(pvad_cols)), ~ rename_to) %>%
    arrange(pvalue) -> res
  
  
  res$cc <- 'NS'
  res[which(abs(res$logFC) > lfcThreshold), 'cc'] <- fc
  res[which(abs(res$padj) <= alpha), 'cc'] <- pv
  res[which(res$padj <= alpha & abs(res$logFC) > lfcThreshold), 'cc'] <- sigfc
  
  up <- 'Up-regulated'
  down <- 'Down-regulated'
  
  res$lfcT <- 'Basal'
  res[which(res$padj < alpha & res$logFC > lfcThreshold), 'lfcT'] <- up
  res[which(res$padj < alpha & res$logFC < -lfcThreshold), 'lfcT'] <- down
  
  res %>%
    # sample_NS(.,1000) %>%
    mutate(cc = factor(cc, levels = c(sigfc, pv, fc, "NS"))) %>%
    mutate(lfcT = factor(lfcT, levels = c(up, down, 'Basal')))
  
}

paste_go <- function(x) { 
  x <- x[!is.na(x)] 
  x <- unique(sort(x))
  x <- paste(x, sep = ';', collapse = ';')
}
