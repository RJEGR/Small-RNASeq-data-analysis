
# 


gene2GO %>% 
  mutate(GOs = strsplit(GOs, ";")) %>%
  unnest(GOs) %>% rename("GO.ID" = "GOs") %>%
  distinct() %>%
  right_join(distinct(DataViz, GO.ID, f1, biotype_best_rank, Term, parentTerm)) %>% 
  count(f1, MajorRNA, Term, biotype_best_rank, parentTerm) %>%
  ggplot(aes(x = MajorRNA, y = parentTerm, fill = n)) + geom_tile(width = 0.5, height = 0.5) + facet_grid(biotype_best_rank ~ f1, scales = "free_x", space = "free_x")
  
  
  

if(!is.null(dev.list())) dev.off()

library(tidyverse)

dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

# DATA <-read_tsv(file.path(dir, "Boostrap_topGO.tsv"))
DATA <-read_tsv(file.path(dir, "Boostrap_topGO_intra_inter.tsv"))
# file.path(dir, "Boostrap_topGO_intra_inter.tsv")

recode_c <- c("CONTRAST_A:-logFC" = "24 hpf:pH 7.6",
  "CONTRAST_B:+logFC" = "110 hpf:pH 8.0",
  "CONTRAST_B:-logFC" = "110 hpf:pH 7.6",
  "CONTRAST_C:+logFC" = "24 hpf:pH 8.0",
  "CONTRAST_C:-logFC"= "110 hpf:pH 8.0",
  "CONTRAST_D:+logFC"= "24 hpf:pH 7.6",
  "CONTRAST_D:-logFC" = "110 hpf:pH 7.6")


DATA <- DATA %>%
  filter(grepl("CONTRAST_C", GROUP)) %>%
  mutate(CONTRAST = dplyr::recode_factor(GROUP, !!!recode_c)) %>%
  mutate(Top = ifelse(!Top %in% c(20, 50, 75, 100), "All", Top)) %>%
  mutate(Top = factor(as.character(Top), levels = c(20, 50, 75, 100, "All"))) %>%
  group_by(CONTRAST, Top) %>%
  mutate(Ratio = size / max(size)) %>%
  filter(Ratio > 0) %>%
  separate(CONTRAST, into = c("f1", "f2"), sep = ":")

DATA %>%
  mutate(classicKS = as.double(classicKS), 
    elimKS = as.double(elimKS)) %>%
  ggplot() + geom_histogram(aes(elimKS))

DATA %>% ungroup() %>% count(f1, f2)

DATA %>%
  mutate(classicKS = as.double(classicKS), 
    elimKS = as.double(elimKS)) %>%
  filter(classicKS <0.05) %>%
  ggplot(aes(x = Top, y = term, fill = Ratio)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_grid(f1 ~f2, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # guides(fill = "none") +
  labs(y = "Biological process", x = "Top Enrichment") +
  scale_fill_viridis_c("Size", option = "inferno", direction = -1)


p <- DATA %>%  
  # filter(as.numeric(classicKS) < 0.05) %>%
  group_by(parentTerm, f1, f2) %>% 
  # count()
  summarise(Ratio = size/sum(size)) %>% 
  ggplot(aes(x = f2, y = parentTerm, fill = Ratio)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  facet_grid(~ f1, switch = "y", scales = "free_y") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  # guides(fill = "none") +
  labs(y = "Biological process", x = "") +
  scale_fill_viridis_c("Enrichment ratio", option = "inferno", direction = -1)

p +
  theme(legend.position = "top",
    # axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    strip.background.x = element_rect(fill = 'grey89', color = 'white'),
    strip.text.y.left = element_text(angle = 0, size = 10, hjust = 1),
    strip.background.y = element_rect(fill = 'white', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 10),
    axis.text.x = element_text(angle = 90, size = 10)) 


##

DataViz <- DATA %>%  
  filter(as.numeric(elimKS) < 0.05) %>%
  filter(Top == 50)

# DataViz %>% count(f1, f2, biotype_best_rank,parentTerm) %>% view()
  
Dataviz <- DataViz %>%
  # mutate(col = -log10(p.adj.ks)) %>%
  group_by(f1, f2, biotype_best_rank) %>%
  mutate(frac = size / max(size)) %>%
  arrange(desc(frac), .by_group = T) %>% 
  mutate(Label = term, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  ungroup() %>% filter(frac > 0.05)

col_recode <- structure(c("#FFC107", "#2196F3","red"), 
  names = c("Intergenic", "Intragenic", "Other ncRNA"))


Dataviz %>%
  # mutate(Label = stringr::str_to_sentence(Label)) %>%
  mutate(f1 = factor(f1, levels = c("24 hpf","110 hpf"))) %>%
  ggplot(aes(y = Label, x = size, color = biotype_best_rank)) + # 
  # facet_wrap(~f1, scales = "free_y", nrow = 2) +
  facet_grid(f1~., scales = "free_y", space = "free_y", switch = "y") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3, alpha = 1) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x), position = "right") +
  labs(y = "microRNA features (Biological process)", x = "Enrichment ratio") +
  theme_bw(base_family = "GillSans", base_size = 12) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1),
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    axis.text.y = element_text(size = 5),
    # axis.ticks.y = element_blank()
    # axis.text.x = element_text(angle = 90, size = 7)
  ) -> p2


p2 <- p2 +
  geom_text(aes(label=parentTerm), x = 0, hjust=0, vjust = 0.1, size = 1, family = "GillSans", color = "gray40")

ggsave(p2, filename = 'TopGO_by_Contrast_3__.png', path = dir, 
  width = 3.7, height = 4, device = png, dpi = 300)


Dataviz %>%
  ungroup() %>%
  count(term,parentTerm, sort = T)

data_text <- Dataviz %>% 
  mutate(Text = gsub("__.+$", "", Label))
  # mutate(Text = stringr::str_to_sentence(Text))


p3 <- p2 + 
  theme( axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  geom_text(data = data_text, 
    aes(label=Text), x = 0, hjust=0, vjust = 0.1, size = 2, family = "GillSans", color = "gray40") 

ggsave(p3, filename = 'TopGO_by_Contrast_3.png', path = dir, 
  width = 3.5, height = 7.5, device = png, dpi = 300)



DATA %>%  
  filter(as.numeric(elimKS) < 0.05) %>%
  filter(Top == 50) %>%
  group_by(f1, f2, biotype_best_rank, parentTerm) %>%
  # mutate(frac = size / max(size)) %>%
  # group_by(f1, f2, biotype_best_rank, parentTerm) %>%
  summarise(size = sum(size)) %>%
  mutate(frac = size / max(size)) %>%
  arrange(desc(frac), .by_group = T) %>% 
  group_by(f1, f2) %>%
  mutate(Label = parentTerm, row_number = row_number(Label)) %>% 
  mutate(Label = factor(paste(Label, row_number, sep = "__"), 
    levels = rev(paste(Label, row_number, sep = "__")))) %>%
  ungroup() %>% 
  # filter(frac > 0.01) %>%
  mutate(f1 = factor(f1, levels = c("24 hpf","110 hpf"))) %>%
  ggplot(aes(y = Label, x = frac, color = biotype_best_rank)) + # 
  # facet_wrap(~f1, scales = "free_y", nrow = 2) +
  facet_grid(f1~., scales = "free_y", space = "free_y", switch = "y") +
  geom_segment(aes(xend = 0, yend = Label), linewidth = 1.3) +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x), position = "right") +
  labs(y = "microRNA features (Biological process)", x = "Enrichment ratio") +
  theme_bw(base_family = "GillSans", base_size = 14) +
  scale_color_manual("",values = col_recode) +
  scale_fill_manual("",values = col_recode) +
  theme(legend.position = "top",
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1),
    # panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    axis.text.y = element_text(size = 5),
    # axis.ticks.y = element_blank()
    # axis.text.x = element_text(angle = 90, size = 7)
  ) -> p4


ggsave(p4, filename = 'TopGO_by_Contrast_4.png', path = dir, 
  width = 5, height = 4, device = png, dpi = 300)

# Back miRNA ID

data_text
