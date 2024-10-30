# 


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

rm(list = ls());

if(!is.null(dev.list())) dev.off()

path <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

f <- file.path(path, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds")

datExpr <- read_rds(f)


PCA <- function(datExpr) {
  
  data <- DESeq2::varianceStabilizingTransformation(round(datExpr))
  
  PCA <- prcomp(t(data), center = T, scale. = FALSE)
  
  PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], Method = basename(f))
  
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  # percentVar <- round(PCA$sdev/sum(PCA$sdev)*100,1)
  
  PCAvar <- data.frame(
    Eigenvalues = PCA$sdev,
    percentVar = percentVar,
    Varcum = cumsum(percentVar),
    Method = basename(f))
  
  return(list(PCAdf, PCAvar))
}

PCAdf <- PCA(datExpr)

Manifest <- data.frame(LIBRARY_ID = colnames(datExpr)) %>%
  mutate(HPF = ifelse(grepl("110", LIBRARY_ID), "100 hpf", "24 hpf")) %>%
  mutate(pH = ifelse(grepl("76", LIBRARY_ID), "pH 7.6", "pH 8.0"))

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
  mutate(LIBRARY_ID = gsub("HR","", LIBRARY_ID)) %>% 
  ggplot(., aes(PC1, PC2)) +
  theme_bw(base_size = 14, base_family = "GillSans") +
  # xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  # ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  xlab(percentVar1) + ylab(percentVar2) +
  # scale_color_grey("") +
  # facet_grid(~ Method) +
  ylim(-20,20) +
  xlim(-20,20) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  # ggrepel::geom_text_repel( family = "GillSans", mapping = aes(label = HPF), size = 3.5, color = "black") +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 3.5, color = "black") +
  theme( legend.position = 'top',
    legend.key.width = unit(0.2, "mm"),
    legend.key.height = unit(0.2, "mm")) -> p

p <- p + 
  ggforce::geom_mark_circle(aes(group = as.factor(HPF), label = HPF), 
    fill = 'grey89', color = NA, label.fontsize = 8, 
    con.cap = unit(0.5, "mm"), #expand = unit(0, "mm"),
    con.size = 0.35, con.type = "straight",
    label.family = "GillSans") +
  geom_point(aes(color = pH), size = 3.7, alpha = 0.7) +
  scale_color_manual("", values = pHpalette)

p <- p + theme(
  panel.grid = element_blank())

ggsave(p, filename = 'PCA.png', path = path, width = 3.5, height = 3.5, device = png, dpi = 500)


PCAdf[[2]] %>%
  # mutate(Method = dplyr::recode_factor(Method, !!!recode_to)) %>%
  group_by(Method) %>%
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
