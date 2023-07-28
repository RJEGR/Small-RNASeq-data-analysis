# 


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
