# Open strucviz results (ViennaRNA format)
# condense the position distribution for each mir 
# write condense distr w/ E-score ex. ..))))..... (-30.70)
# visualize 
# example, using forna (https://github.com/ViennaRNA/forna)
 # which use force-directed graph layout 


rm(list = ls());

if(!is.null(dev.list())) dev.off()

dir <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

print(.KnownRNAsDB <- read_rds(paste0(dir, "RNA_LOCATION_MIR_DB.rds")))

KnownRNAsDB <- .KnownRNAsDB %>%
  mutate(KnownRNAs = strsplit(KnownRNAs, ";")) %>%
  dplyr::select(Name, KnownRNAs, MajorRNA) %>%
  unnest(KnownRNAs) %>%
  mutate(KnownRNAs = gsub("_[3-5]p$","", KnownRNAs)) %>%
  mutate(KnownRNAs = ifelse(is.na(KnownRNAs), "Novel (47)", "Known (70)")) %>%
  distinct(Name, KnownRNAs) 

dir <- "/Users/cigom/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/Outputs/strucVis"

f <- list.files(path = dir, pattern = "txt", full.names = T)

# f <- f[1]


library(tidyverse)

# read only energy of folding
read_escore <- function(f) {
  
  head(x <- read_lines(f, n_max =1, skip = 6))
  
  # this input for force-directed graph layout 
  # split_str <- sapply(strsplit(x, "-"), `[`, 1)
  # http://rna.tbi.univie.ac.at/forna/
  
  split_str <- sapply(strsplit(x, " "), `[`,2)

  escore <- as.double( gsub("\\(|\\)", "", split_str))
  
  Name <- gsub(".txt", "", basename(f))
  
  out <- data.frame(  Name, escore)
  
  return(out)
}

# read_escore(f[1])

read_struc <- function(f) {
  
  split_str <- read_lines(f, skip = 7)
  
  seqs <- sapply(strsplit(split_str, "len:"), `[`, 1)
  
  dim(seqs <- stringr::str_split(seqs, pattern = "", simplify = T))
  
  split_str <- sapply(strsplit(split_str, "len:"), `[`, 2)
  
  len <- sapply(strsplit(split_str, "al:"), `[`, 1)
  
  al <- sapply(strsplit(split_str, "al:"), `[`, 2)
  
  len <- as.numeric(len)
  
  al <- as.numeric(al)
  
  # out <- data.frame(seqs, len, al) %>% as_tibble()
  
  # return(out)
  
  # or as single vector
  
  na_fill <- function(x) { gsub("[.]|[[:space:]]", NA, x)} 
  
  nuc_freq <- function(x) {sum(!is.na(na_fill(x)))}
  
  out <- data.frame(seqs)
  
  out <- apply(out, 2, nuc_freq)
  
  Precs <- max(al)/sum(al)
  
  Name <- gsub(".txt", "", basename(f))
  
  out <- data.frame(Name = Name,
    Reads = sum(al), precision =  max(al)/sum(al), n_seqs = nrow(seqs),
    Freq = out) %>%
    rowid_to_column(var = "pos") 
  
  return(out)
  
}

escoresdf <- lapply(f, read_escore)
escoresdf <- do.call(rbind, escoresdf)

pdens <- escoresdf %>%
  left_join(KnownRNAsDB) %>%
  ggplot(aes(y = KnownRNAs, x = escore, fill = stat(x))) +
  # facet_grid(x ~ ., scales = "free_y", space = "free", switch = "y") +
  facet_wrap(~ KnownRNAs, nrow = 2, scales = "free_y") +
  labs(y = "", x = "Minimum free-energy (MFE)") +
  scale_fill_viridis_c(option = "C") +
  # xlim(0,1) +
  ggridges::geom_density_ridges_gradient(
    jittered_points = T,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 0.5, alpha = 0.2) +
  # scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_classic(base_family = "GillSans", base_size = 8) +
  theme(
    strip.text = element_text(color = "black",hjust = 1),
    strip.background = element_blank(),
    legend.position = 'none',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 8),
    # axis.line.x = element_blank(),
    panel.grid.major = element_blank())


out <- lapply(f, read_struc)

out <- do.call(rbind, out) %>% as_tibble()

out %>% distinct(Name, Reads, precision, n_seqs)

out <- out %>% mutate(col = Freq/n_seqs) %>%  left_join(KnownRNAsDB) 

# out %>% distinct(n_seqs, KnownRNAs) %>% group_by(KnownRNAs) %>% summarise(sum(n_seqs))

lo = floor(min(out$col))
up = ceiling(max(out$col))
mid = (lo + up)/2

p <- out %>% 
  mutate(Freq = ifelse(Freq > 0, Freq, NA)) %>%
  drop_na(Freq) %>%
  ggplot(aes(x = pos, y = Name, color = col)) + 
  # geom_tile(size = 0.5) +
  geom_point(shape = 15, size = 0.25) +
  facet_grid(KnownRNAs ~ ., scales = "free") +
  scale_color_gradient2(
    low = "#FFFC00", high = "red", mid = "#6C379E", 
    # low = "blue", high = "red", mid = "white", 
    na.value = "white", midpoint = mid, limit = c(lo, up),
    name = NULL) +
  scale_x_continuous(position = "top", limits = c(0,100)) +
  labs(y = "microRNA precursor", x = "") +
  theme_bw(base_family = "GillSans", base_size = 10) +
  guides(color = guide_colorbar(
    # barwidth = unit(2, "in"),
    # barheight = unit(0.1, "in"), 
    title = "% Depth of Coverage",
    title.position  = "top",
    title.theme = element_text(size = 8, family = "GillSans", hjust = 1),
    label.position = "bottom",
    alignd = 0.5,
    ticks.colour = "black", ticks.linewidth = 0.35,
    frame.colour = "black", frame.linewidth = 0.35,
    label.theme = element_text(size = 8, family = "GillSans"))) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1),
    legend.position = 'bottom',
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.12, "cm"),
    # panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "inside",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p

path_out <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

library(patchwork)

ps <- p +  plot_spacer() + pdens + plot_layout(widths = c(0.4, -0.045, 0.25))

# ps

ggsave(ps, filename = 'FIGURE_1_PANEL_SUP1.png', path = path_out, 
  width = 3.5, height = 5, device = png, dpi = 500)


top <- out %>% 
  filter(Freq > 0) %>%
  group_by(pos, KnownRNAs) %>%
  summarise(Freq = sum(Freq)) %>%
  mutate(facet = "Sequence length") %>%
  ggplot(aes(x = pos, y = Freq, color = KnownRNAs)) + 
  facet_grid(~ facet) +
  labs(x = "Sequence length", y = "", color = "") +
    # geom_point(size = 0.5) +
  geom_step() +
  scale_x_continuous(position = "bottom", limits = c(0,100)) +
  scale_color_grey() +
  theme_classic(base_family = "GillSans", base_size = 10) +
  theme(
    strip.background = element_rect(fill = 'white', color = 'white'),
    strip.text = element_text(color = "black",hjust = 1),
    legend.position = 'top',
    legend.key.width = unit(0.2, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank())



ps <- top /  plot_spacer() / p + plot_layout(heights = c(0.07, -0.055, 0.5))


ggsave(ps, filename = 'FIGURE_1_PANEL_SUP.png', path = path_out, 
  width = 2.5, height = 5.7, device = png, dpi = 500)


# ps <- ps +  plot_spacer() + pdens + plot_layout(widths = c(0.4, -0.045, 0.25))




# 
# apply(m, 2, nuc_freq) %>%
#   as_tibble()

cols <- m  %>% select_if(is.character) %>% colnames()

# which the majorRNA

m %>% 
  rowid_to_column(var = "seq") %>%
  pivot_longer(all_of(cols), names_to = "pos", values_to = "nuc") %>%
  mutate(pos = gsub("^X", "", pos), pos = as.numeric(pos)) %>%
  mutate(nuc = na_fill(nuc)) %>%   drop_na(nuc) %>%
  mutate(nuc = stringr::str_to_upper(nuc)) %>%
  ggplot(aes(y = seq, x = pos, fill = nuc, alpha = al)) + geom_tile()
# Freq of nuc.

m %>% 
  rowid_to_column(var = "seq") %>%
  pivot_longer(all_of(cols), names_to = "pos", values_to = "nuc") %>%
  mutate(pos = gsub("^X", "", pos), pos = as.numeric(pos)) %>%
  mutate(nuc = na_fill(nuc)) %>% 
  drop_na(nuc) %>%
  count(pos, nuc) %>%
  mutate(nuc = stringr::str_to_upper(nuc)) %>%
  ggplot(aes(x = pos, y = n, fill = nuc)) + geom_col() +
  xlim(1,100)

# alignment dist

m %>% 
  rowid_to_column(var = "seq") %>%
  pivot_longer(all_of(cols), names_to = "pos", values_to = "nuc") %>%
  mutate(pos = gsub("^X", "", pos), pos = as.numeric(pos)) %>%
  mutate(nuc = na_fill(nuc)) %>% 
  drop_na(nuc) %>%
  group_by(pos) %>%
  summarise(n = nuc_freq(nuc)) %>%
  # group_by(seq,pos) %>% summarise(al = sum(al), n = n()) %>% view()
  # mutate(al = al*1E-6) %>%
  ggplot(aes(x = pos, y = n)) + 
  geom_point(size = 0.5) +
  geom_step() +
  xlim(1,100)

