# AFTER PREPARING DATA USING respirometries.R
# LOAD DATA
# GENERATE DATAVIZ
# RICARDO GOMEZ-REYES, 2023

# Note: values were converted to pmol by: 
# From (μmol O2⋅larva− 1 ⋅ h− 1 L-1)
# To (μmol O2⋅larva− 1 ⋅ h− 1).
# Then ((pmol O2⋅larva− 1 ⋅ h− 1))

# Vol of chamber: 1700 μL OR 0.0017 L

# df_stats <- df_stats %>% 
#   mutate(r_ind_adj = r_ind_adj * 0.0017) %>%
#   mutate(r_ind_adj = r_ind_adj * 1E6)



rm(list = ls());

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE)

library(tidyverse)
library(rstatix)

path_out <- '~/Documents/MIRNA_HALIOTIS/'

.x <- read_rds(paste0(path_out, 'resp_rates.rds'))

df_stats <- .x[[1]]

.stats <- .x[[2]]

# pHpalette <- c(`7.6`="#d73027", `7.8`= "#abdda4",`8.0`= "#4575b4", `8`= "#4575b4")

pHpalette <- c(`7.6`="#ad1f1f", `7.8`= "#abdda4",`8.0`= "#4575b4", `8`= "#4575b4")

pHpalette <- c(`7.6`="black", `7.8`= "#abdda4",`8.0`= "gray68", `8`= "gray68")



recode_to <- structure(c("pH 8.0", "pH 7.6"), names = c("8", "7.6"))

recode_hpf <- structure(c("B) 24 hpf", "C) 110 hpf"), names = c("24", "108"))

# View

df_stats %>% select(r_ind_adj) %>%
  # filter(pH != "7.8") %>%
  dplyr::mutate(hpf = dplyr::recode(hpf, !!!c("108" = "110"))) %>%
  rstatix::get_summary_stats(type = 'mean_sd') # %>% view()


# 1) (REPLACED W/ DONUTS)

ylabs <- expression("RM (pmol"~O[2]~Ind^-1~h^-1*")")

#ylabs <- expression("Tasa de resp. (pmol"~O[2]~Ind^-1~h^-1*")")

lim <- c(-10, 140)

breaks <- seq(0, lim[2], by = 20)

df_stats %>% select(r_ind_adj) %>%
  filter(pH != "7.8") %>%
  dplyr::mutate(hpf = dplyr::recode(hpf, !!!c("108" = "110"))) %>%
  rstatix::get_summary_stats(type = 'mean_sd') %>% 
  mutate(facet = "A) Larval development") %>%
  # mutate(facet = "A) Desarrollo larval") %>%
  mutate(ymin = mean-sd, ymax = mean+sd) %>%
  ggplot(aes(x = hpf, y = mean, color = pH, group = pH)) +
  facet_grid( ~ facet ) +
  geom_point(position = position_dodge(width = 0), size = 3, alpha = 0.5) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
    width = 0.1, position = position_dodge(width = 0)) +
  geom_path(position = position_dodge(width = 0), linewidth = 1) +
  scale_color_manual("", values = pHpalette) +
  labs(y = ylabs, x = "Time (hpf)") +
  scale_y_continuous(breaks = breaks, limits = lim) -> pleft1



pleft1 + theme_bw(base_family = "GillSans", base_size = 14) + 
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none") -> pleft1

# ggsave(pleft,
#   filename = 'oxygen_rate_facet_by_hpf.png', path = ggsavepath,
#   width = 3, height = 2.5)


# 1) 
recode_to_pH <- structure(c("A) pH 8.0", "B) pH 7.8", "C) pH 7.6"), names = c("8","7.8","7.6"))


donut_df <- df_stats %>% select(r_ind_adj) %>%
  # filter(pH != "7.8") %>%
  dplyr::mutate(hpf = dplyr::recode(hpf, !!!c("108" = "110"))) %>%
  rstatix::get_summary_stats(type = 'mean_sd')
  

donut_df <- donut_df %>%
  group_by(pH) %>%
  mutate(Frac = mean/sum(mean)) %>%
  dplyr::mutate(pH = dplyr::recode(pH, !!!recode_to_pH))

donut_df %>% summarise(Frac = sum(Frac))

hsize <- 1.5

pleft <- donut_df %>%
  mutate(label = paste0(round(Frac*100)," %")) %>%
  mutate(x = hsize) %>%
  ggplot(aes(x = x, y = Frac, fill = hpf)) +
  facet_grid(~ pH) +
  geom_col(color = "white") + # aes(alpha = Frac)
  geom_text(aes(label = label),
    position = position_stack(vjust = 0.5), family = "GillSans") +
  coord_polar(theta = "y") +
  scale_fill_brewer("HPF", palette = "GnBu") +
  theme_bw(base_family = "GillSans", base_size = 14) +
  xlim(c(0.2, hsize + 0.5)) +
  theme(
    legend.position = "top",
    # panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank())

# 1.2) barplot
# SUMMARISE AVERAGE BY REPLICATE

df_stats %>% ungroup() %>% count(pH, Spot)

# 4 spots == 4 stages of development (24, 48, 60 and 110)

DF <- df_stats %>% 
  # filter(Spot %in% which_spot) %>%
  group_by(pH, Spot) %>%
  summarise(y = sum(r_ind_adj), n = n()) %>%
  filter(n == 4)

# BECAUSE PARAMETRIC AND HOMOCELASTICITY

# IS DIFF BETWEEN STAGES?

df_stats %>% 
  # group_by(pH) %>%
  ungroup() %>%
  rstatix::anova_test(r_ind_adj ~ hpf) %>%
  add_significance()

# PRIOR
# TEST IF DIFFERENCE B/ TOTAL AMOUNT OF O CONSUMED:

df_stats %>% 
  group_by() %>%
  rstatix::anova_test(r_ind_adj ~ pH) %>%
  # rstatix::kruskal_test(r_ind_adj ~ pH) %>% # if non param.
  add_significance()

# UNNECESARY, BUT POSTERIOR

df_stats %>%
  ungroup() %>%
  # tukey_hsd(r_ind_adj ~ pH) %>%
  pairwise_t_test(r_ind_adj ~ pH, ref.group = '8') %>%
  # pairwise_wilcox_test(r_ind_adj ~ pH,  conf.level = 0.95) %>%
  adjust_pvalue() %>%
  add_significance()

df_stats %>%
  ungroup() %>%
  # group_by(hpf) %>%
  # tukey_hsd(r_ind_adj ~ pH) %>%
  pairwise_t_test(r_ind_adj ~ pH) %>%
  # pairwise_wilcox_test(r_ind_adj ~ pH,  conf.level = 0.95) %>%
  adjust_pvalue() %>%
  add_significance()

recode_to_x <- structure(c("pH 8.0", "pH 7.8", "pH 7.6"), names = c("8","7.8","7.6"))
  
# ylabs <- expression("Metabolic rate (pmol"~O[2]~Ind^-1~h^-1*")")

pbar <- DF %>%
  summarise(mean = mean(y), sd = sd(y), n = n()) %>%
  dplyr::mutate(pH = dplyr::recode(pH, !!!recode_to_x)) %>%
  mutate(ymin = mean-sd, ymax = mean+sd, p.adj.signif = c("", "ns", "ns")) %>%
  mutate(facet = "D) Total oxygen consumed") %>%
  # mutate(facet = "D) Consumo total de O2") %>%
  ggplot(aes(x = pH, y = mean,  ymin = ymin, ymax = ymax)) +
  facet_grid(~ facet) +
  geom_col(position = position_dodge(0.6), color = 'black', fill = 'grey89') +
  geom_errorbar(width = 0.15, position = position_dodge(0.6)) +
  geom_text(aes(y = ymax + 2, label= p.adj.signif), 
      size = 3.5, family = 'GillSans',fontface = "bold",
      vjust= 0, color="black", position=position_dodge(width = .6)) +
  # scale_y_continuous(breaks = seq(0, 0.11, by = 0.02), limits = c(0,0.11)) +
  theme_bw(base_family = "GillSans", base_size = 14) + 
  theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none") +
  labs(y = ylabs, x = "") 

ggsave(pbar, filename = 'RESPIRATION_RATES_FACET_3.png', 
  path = path_out, width = 3, height = 3, device = png, dpi = 300)



# DF %>%
#   dplyr::mutate(pH = dplyr::recode(pH, !!!recode_to_x)) %>%
#   ggplot(aes(x = pH, y = y)) +
#   geom_col(position = position_dodge2(0.6), color = 'black', fill = 'white') +
#   # geom_errorbar(width = 0.15, position = position_dodge(0.6)) +
#   theme_bw(base_family = "GillSans", base_size = 14) + 
#   theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
#     panel.border = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     legend.position = "none")

# 2)

fun.data.trend <- "mean_sdl" # "mean_cl_boot", "mean_se"
  
plotdf <- df_stats %>% 
  filter(pH %in% names(recode_to)) %>%
  filter(hpf %in% names(recode_hpf)) %>%
  dplyr::mutate(hpf = dplyr::recode(hpf, !!!recode_hpf)) %>%
  dplyr::mutate(pH = dplyr::recode(pH, !!!recode_to)) 

plotdf %>%
  ungroup() %>%
  rstatix::anova_test(r_ind_adj ~ hpf) %>%
  # rstatix::kruskal_test(r_ind_adj ~ pH) %>% # if non param.
  add_significance()
  

stats <- .stats %>% 
  filter(group1 == "8" & group2 == "7.6") %>%
  filter(hpf %in% names(recode_hpf)) %>%
  dplyr::mutate(hpf = dplyr::recode(hpf, !!!recode_hpf)) %>%
  dplyr::mutate(group1 = dplyr::recode(group1, !!!recode_to)) %>%
  dplyr::mutate(group2 = dplyr::recode(group2, !!!recode_to)) %>%
  mutate(xmin = "pH 8.0", xmax = "pH 7.6", )
  # select(hpf, group2, p.adj.signif) %>% dplyr::rename("pH" = "group2") %>%
  # right_join(df_stats)


plotdf %>%
  ggplot(aes(x = pH, y = r_ind_adj, group = hpf)) +
  facet_wrap(~ hpf) +
  geom_jitter(aes(color = pH, fill = pH), position=position_jitter(w=0.1,h=0), size = 1, alpha = 1) +
  stat_summary(fun.data = fun.data.trend, colour = "black", linewidth = 0.7, size = 0.5, alpha = 0.7) +
  stat_summary(fun = mean, geom = "line", colour = "black", alpha = 0.7) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  labs(y = ylabs, x = "Assay") +
  scale_color_manual("", values = pHpalette) + # c("#4575b4", "#d73027")
  scale_fill_manual("", values = pHpalette) + # c("#4575b4", "#d73027")
  # scale_y_continuous(breaks = seq(0, 0.11, by = 0.02), limits = c(0,0.11))
  scale_y_continuous("",breaks = breaks, limits = lim) -> pright


pright <- pright + theme(strip.background = element_rect(fill = 'grey89', color = 'white'),
  panel.border = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  legend.position = "none") 


pright <- pright + 
  ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
    remove.bracket = F, tip.length = 0.01,  hide.ns = T, family = "GillSans") 

# ps <- ps + geom_text(aes(y = ymax + 0.04, label= p.adj.signif), 
#   size = 6, family = 'GillSans',fontface = "bold",
#   vjust= 0, color="black", position=position_dodge(width = .7))


library(patchwork)


# p <- pleft + plot_spacer() + pright + patchwork::plot_layout(widths = c(0.8,-0.15, 1.2))

# p <- pleft / plot_spacer() / pright + patchwork::plot_layout(heights = c(1,-0.15, 1))

# p <- pleft / plot_spacer() / pright + patchwork::plot_layout(widths = c(3, -1.2, 3), heights = c(20,0,4))

cowplot::plot_grid(pleft, pbar, ncol=2, align='hv', nrow = 1, rel_widths = c(1, .6), rel_heights = c(1.25, 0.2))

# wrap_plots(pleft, pbar)

ggsave(pleft, filename = 'RESPIRATION_RATES_FACET_2.png', 
  path = path_out, width = 5, height = 3, device = png, dpi = 300)

ggsave(pbar, filename = 'RESPIRATION_RATES_FACET_3.png', 
  path = path_out, width = 3, height = 3, device = png, dpi = 300)

ggsave(pright, filename = 'RESPIRATION_RATES_FACET_4.png', 
  path = path_out, width = 5, height = 3, device = png, dpi = 300)


pleft1 | pright

P <- pleft1 + plot_spacer() + pright + patchwork::plot_layout(widths = c(0.8,-0.15, 1.2))

ggsave(P, filename = 'RESPIRATION_RATES_FACET_EN.png', 
  path = path_out, width = 6, height = 3, device = png, dpi = 300)
