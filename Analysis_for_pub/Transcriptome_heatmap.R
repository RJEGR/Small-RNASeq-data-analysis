# Generate heatmap for suppl material
# using transcriptome data
# Also try calculate transcriptional noise from lm


rm(list = ls());

if(!is.null(dev.list())) dev.off()


library(tidyverse)

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/FUNCTIONAL_MIR_ANNOT/"

print(TARGETDB <- read_rds(paste0(dir,"SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds")))

Mantle_sam <- c("SRR8956768", "SRR8956769")

# remove mantle samples

.TARGETDB <- TARGETDB %>% dplyr::select(!starts_with(Mantle_sam)) 

is_na <- function(x) ifelse(is.na(x), 0, x)

keep_expressed <- .TARGETDB %>% 
  dplyr::select(starts_with("SRR")) %>% 
  mutate(across(where(is.double), ~is_na(.))) %>% 
  rowSums()

sum(keep_expressed <- keep_expressed > 1) # 131

.TARGETDB <- .TARGETDB[keep_expressed,]

which_cols <- .TARGETDB %>% select_if(is.double) %>% colnames()

nrow(.TARGETDB <- .TARGETDB %>%
  group_by(gene_id, description) %>%
  summarise_at(vars(all_of(which_cols)), sum))


COUNT <- .TARGETDB %>%
  ungroup() %>%
  select_if(is.double) %>% as("matrix")

rownames(COUNT) <- .TARGETDB$gene_id

LIBRARY_ID <- colnames(COUNT)


dir <- "~/Documents/MIRNA_HALIOTIS/"

wd <- paste0(dir, "RNA_SEQ_ANALYSIS/ASSEMBLY/STRINGTIE/REFBASED_MODE/")

Manifest <- read_csv(list.files(path = wd, pattern = "METADATA_RNASEQ.csv", full.names = T))

Manifest[is.na(Manifest)] <- "Low CO2"

Manifest <- Manifest[!grepl("Mantle", Manifest$Isolated),]

Manifest <- Manifest %>% arrange(match(LIBRARY_ID, colnames(count)))

any(colnames(COUNT) == Manifest$LIBRARY_ID) # sanity check

Manifest <- mutate_if(Manifest, is.character, as.factor)


tmp.dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

data.frame(COUNT) %>% 
  as_tibble(rownames = 'gene_id') %>%
  pivot_longer(cols = colnames(COUNT), values_to = 'Read_expression', names_to = "LIBRARY_ID") %>%
  left_join(Manifest, by = "LIBRARY_ID") %>% filter(Read_expression > 0) %>%
  write_tsv(file.path(tmp.dir, "Transcriptome_read_expression.tsv"))


# PCA


PCA <- function(datExpr) {
  
  # require(DESeq2)
  
  data <- DESeq2::varianceStabilizingTransformation(round(datExpr))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar)
    # Method = basename(f)
    )
  
  return(list(PCAdf, PCAvar))
}

PCAdf <- PCA(COUNT)


percentVar1 <- paste0(PCAdf[[2]]$percentVar[1], "% ,",PCAdf[[1]]$percentVar[1])
percentVar1 <- paste0("PC1, VarExp: ", percentVar1, "%")


percentVar2 <- paste0(PCAdf[[2]]$percentVar[2], "% ,",PCAdf[[1]]$percentVar[2])
percentVar2 <- paste0("PC2, VarExp: ", percentVar2, "%")


PCAdf_ <- PCAdf[[1]]

# recode_to <- c("Trinity_counts.txt", "Rnaspades_counts.txt")
# recode_to <- structure(c("Trinity", "Rnaspades"), names = recode_to)

pHpalette <- c(`pH 7.6`="black",`pH 8.0`= "gray78", `pH 8`= "gray78")


PCAdf_ %>%
  # mutate(Method = dplyr::recode_factor(Method, !!!recode_to)) %>%
  mutate(LIBRARY_ID = rownames(.)) %>% 
  left_join(Manifest, by = "LIBRARY_ID") %>%
  # mutate(LIBRARY_ID = gsub("HR","", LIBRARY_ID)) %>% 
  ggplot(., aes(PC1, PC2)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  # xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  # ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  xlab(percentVar1) + ylab(percentVar2) +
  # scale_color_grey("") +
  # facet_grid(~ Method) +
  # ylim(-20,20) +
  # xlim(-20,20) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = HPF), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme( legend.position = 'top',
    legend.key.width = unit(0.2, "mm"),
    legend.key.height = unit(0.2, "mm")) -> p

p + geom_point(aes(color = CONTRAST_A)) + 
  geom_text(aes(label = CONTRAST_B)) +
  theme(
    panel.grid = element_blank())

p <- p + 
  # ggforce::geom_mark_circle(aes(group = as.factor(dpf), label = dpf), 
  #   fill = 'grey89', color = NA, label.fontsize = 8, 
  #   con.cap = unit(0.5, "mm"), #expand = unit(0, "mm"),
  #   con.size = 0.35, con.type = "straight",
  #   label.family = "GillSans") +
  geom_point(aes(color = CONTRAST_B), size = 3.7, alpha = 0.7) +
  scale_color_manual("", values = pHpalette)


PCAdf[[2]] %>%
  # mutate(Method = dplyr::recode_factor(Method, !!!recode_to)) %>%
  # group_by(Method) %>%
  mutate(Dim = row_number()) %>%
  ggplot(., aes(y = percentVar, x = as.factor(Dim))) +
  geom_col(position = position_dodge2(width = 0.5), width = 0.5, fill = "black") + # fill = "gray78", color = "gray78"
  geom_line(aes(y = Varcum), group = 1) +
  geom_point(aes(y = Varcum), size = 1.5) +
  # labs(x = "Component Number", y = "Eigenvalue") +
  labs(x = "Principal component", y = "Fraction variance explained (%)") +
  scale_fill_grey("") +
  scale_color_grey("") +
  guides(color=guide_legend(nrow = 1)) +
  theme_bw(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top", 
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()) 


# ggsave(p, filename = 'PCA.png', path = path, width = 3.5, height = 3.5, device = png, dpi = 500)

# heatmap ----

z_scores <- function(x) {(x-mean(x))/sd(x)}

# Noise_score <- function(x) {1 - x/sum(x)}

barplot(colSums(COUNT))

M <- DESeq2::varianceStabilizingTransformation(round(COUNT))

# barplot(colSums(M))

M <- apply(M, 2, z_scores)

# M <- apply(M, 2, Noise_score);h <- heatmap(M, keep.dendro = T)

h <- heatmap(t(M), keep.dendro = T)

hc_samples <- as.hclust(h$Colv)
plot(hc_samples)
hc_order <- hc_samples$labels[h$colInd]

hc_genes <- as.hclust(h$Rowv)
plot(hc_genes)
order_genes <- hc_genes$labels[h$rowInd]


# DATA <- data.frame(M) %>% 
#   as_tibble(rownames = 'LIBRARY_ID') %>%
#   pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "gene_id") %>%
#   left_join(Manifest, by = "LIBRARY_ID") %>%
#   mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))

DATA <- data.frame(M) %>% 
  as_tibble(rownames = 'gene_id') %>%
  pivot_longer(cols = colnames(M), values_to = 'fill', names_to = "LIBRARY_ID") %>%
  left_join(Manifest, by = "LIBRARY_ID") %>%
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order)))


labels <- .TARGETDB %>% select(gene_id, description)
# arrange(desc(Name)) %>%
# group_by(Name) %>%
# mutate(Label = Name, row_number = row_number(Label)) %>%
# mutate(Label = factor(paste(Label, row_number, sep = "__"), 
#   levels = rev(paste(Label, row_number, sep = "__"))))


labels <- structure(labels$description, names = labels$gene_id)

labels <- labels[match(order_genes, names(labels))]

identical(order_genes, names(labels))

order_genes <- labels

order_genes <- gsub(", transcript variant X*", "", order_genes)

lo = floor(min(DATA$fill))
up = ceiling(max(DATA$fill))
mid = (lo + up)/2


library(ggh4x)

DATA %>%
  mutate(hpf = paste0(dpf*24, " hpf")) %>%
  ggplot(aes(x = LIBRARY_ID, y = gene_id, fill = fill)) +  
  geom_tile(color = 'white', linewidth = 0.025, width = 1) +
  # facet_grid(~ hpf+pH, scales = "free_x") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  ggh4x::scale_y_dendrogram(hclust = hc_genes, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = order_genes, label_size = 5, label_family = "GillSans")) +
  theme_bw(base_size = 10, base_family = "GillSans") +
  labs(x = '', y = '') +
  theme(
    legend.position = "top",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    # strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) -> P

P <- P + 
  ggh4x::facet_nested(~ dpf, nest_line = T, scales = "free_x", space = "free_x") +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = 'gray89', color = 'white')) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(2.5, "in"),
    barheight = unit(0.07, "in"), 
    label.position = "top",
    title = "Row z-score",
    title.position  = "top",
    title.theme = element_text(size = 7, family = "GillSans", hjust = 1),
    alignd = 0.8,
    ticks.colour = "black", ticks.linewidth = 0.25,
    frame.colour = "black", frame.linewidth = 0.25,
    label.theme = element_text(size = 7, family = "GillSans")))

P

dir_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

ggsave(P, filename = 'Target_transcriptome_heatmap.png', path = dir_out, 
  width = 8, height = 8, device = png, dpi = 300)

# Prep summary table ----

Noise_score <- function(x) {1 - x/sum(x)}

M <- DESeq2::varianceStabilizingTransformation(round(COUNT))

N <- which(rownames(COUNT) %in% "LOC124135816")
COUNT[N,]
barplot(Noise_score(COUNT[N,]))

heatmap(apply(COUNT, 2, Noise_score))
  
tail(sort(rowSums(COUNT)))

colSums(COUNT)


# Trancriptome noise analysis
# Ex.
# Transcriptional noise was measured as the percent coefficient of variation of gene expression [49] and used as the response variable in the following multiple linear regression model: log10(transcriptional noise) ~ gene body methylation + log2(gene expression) + log10(gene length) + Wolbachia infection status

# transcriptional noise (measured by the coefficient of variation of expression [49]) as the response variable. https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1008397

M <- DESeq2::varianceStabilizingTransformation(round(COUNT))

EXP <- data.frame(Read_exp = log10(rowSums(COUNT)), Noise = z_scores(log10(rowSums(COUNT))))
EXP <- EXP %>% as_tibble(rownames = "gene_id")

cor(EXP$Read_exp, EXP$Noise)

EXP %>% ggplot(aes(Read_exp, Noise)) + geom_point()

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

EXP <- DB %>% 
  mutate(gene_id = gene_id) %>%
  # mutate(GENESET = ifelse(is.na(STRINGID), gene_id, STRINGID)) %>%
  drop_na(gene_id) %>% 
  dplyr::count(gene_id, sort = T) %>%
  arrange(desc(n)) %>%
  right_join(EXP)
  
EXP %>% ggplot(aes(Read_exp, n)) + geom_point()

cor(EXP$Read_exp, EXP$n)

EXP <- DATA %>% 
  filter(fill > 2.5) %>% distinct(dpf, gene_id) %>% 
  right_join(EXP) %>%
  mutate(dpf = ifelse(is.na(dpf), "Basal", as.character(dpf))) %>%
  mutate_if(is.character, as.factor)


EXP %>%
  ggplot(aes(x = gene_id, y = Read_exp)) + 
  facet_grid(dpf ~.) +
  geom_col(aes(fill = n))

# Model
# log10(transcriptional noise) ~ microRNA_degree + log2(microRNA expression) + Trait

formula <- "Read_exp ~ n + Noise + Isolated"

fit_model <- function(data, formula) {
  
  
  require(tidymodels)
  
  
  formula <- formula(formula)
  
  data_rec <- recipe(formula, data = data) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_dummy(all_nominal_predictors())
  
  
  
  model_spec <- linear_reg() %>% set_engine("lm")
  
  data_workflow <- workflow() %>%
    add_recipe(data_rec) %>%
    add_model(model_spec)
  
  # Split training and test model set
  
  set.seed(123)
  
  data_split <- initial_split(data, prop = 0.8)
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  data_fit <- data_workflow %>% fit(data = data_train)
  
  # data_predictions <- data_fit %>%
  #   predict(new_data = data_test) %>%
  #   bind_cols(data_test)
  
  results <- data_fit %>%
    extract_fit_parsnip() %>%
    tidy()
  
  
  return(results)
}

fitdf <- fit_model(EXP, formula)


fitdf <- fitdf %>%
  mutate(star = ifelse(p.value <.001, "***", 
    ifelse(p.value <.01, "**",
      ifelse(p.value <.05, "*", "ns"))))

fitdf %>%
  filter(!grepl("Intercept",term)) %>%
  mutate(term = gsub("Isolated_X[1-9].","", term)) %>%
  # mutate(term = gsub("[.]"," ", term)) %>%
  mutate(y_star = estimate + (0.2+std.error)*sign(estimate)) %>%
  mutate(xmin = estimate-std.error, xmax = estimate+std.error) %>%
  ggplot(aes(y = term, x = estimate)) + # color = Intercept
  # facet_grid(~ Intercept) +
  geom_point(size = 0.5) +
  geom_text(aes(x = y_star, label = star),  
    vjust = 0.75, hjust = 0, size= 3.5, 
    color="black", 
    # position=position_dodge(width = .5),
    family =  "GillSans") +
  geom_errorbar(aes(xmin = xmin, xmax = xmax), 
    width = 0.1, alpha = 0.3, color = "black",
    # position=position_dodge(width = .5)
  ) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color = "black") +
  labs(
    y = "",
    x = "Effect Size") +
  # scale_x_continuous(breaks = seq(-1, 1, by = 0.5), limits = c(-1,1)) +
  theme_classic(base_family = "GillSans", base_size = 12) 
