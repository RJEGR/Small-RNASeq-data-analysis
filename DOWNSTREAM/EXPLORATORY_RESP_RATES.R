
rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

# library(Biostrings)
library(rstatix)
library(tidyverse)

mtd <- read_tsv(list.files(path = '~/Documents/MIRNA_HALIOTIS/CLEAN_INPUT/', pattern = 'METADATA', full.names = T))

path_resp <- '~/Documents/MIRNA_HALIOTIS/'

load(paste0(path_resp, 'resp_rates.Rdata'))

namespH <- c("8", "7.6")

recodeL <- c('Control', 'Low')

level_key <- structure(recodeL, names = namespH)

df_stats <- df_stats %>% filter(pH != 7.8) %>% 
  mutate(pH = droplevels(pH)) %>%
  mutate(pH = recode_factor(pH, !!!level_key))
  # mutate(hpf = recode_factor(hpf, !!!structure("110", names = "108")))

# rstatix::get_summary_stats(df_stats, type = "mean_se")

# 1) Parametric priori test----

df_stats %>% 
  group_by(hpf) %>%
  rstatix::welch_anova_test(r_ind_adj ~ pH) %>%
  # rstatix::anova_test(r_ind_adj ~ pH) %>%
  add_significance() -> prior.stats

# 2) Parametric posteriori test ----

df_stats %>%
  group_by(hpf) %>%
  tukey_hsd(r_ind_adj ~ pH) %>%
  # pairwise_t_test(r_ind_adj ~ pH)
  adjust_pvalue() %>%
  add_significance() -> post.test

post.test %>% add_xy_position(x = "pH") -> stats

title <- get_pwc_label(stats)

subtitle <- get_test_label(prior.stats, detailed = F)

# plor as boxplot
# o2/Ind (Umol h-1 L-1)

ylabs <- expression("Rate ("*O[2]~mu*"mol"~h^-1~L^-1~Ind^-1*")")

df_stats %>% 
  ggplot(aes(x = pH, y = r_ind_adj, group = pH)) +
  facet_grid(~ hpf) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  stat_boxplot(aes(color = pH),
    geom ='errorbar', width = 0.3) +
  geom_boxplot(aes(color = pH, fill = pH), 
    width = 0.3, outlier.alpha = 0) +
  labs(y = ylabs) -> psave


psave + theme(strip.background = element_rect(fill = 'grey', color = 'white'),
  panel.border = element_blank(),
  legend.position = "top") -> psave

psave + 
  ggpubr::stat_pvalue_manual(stats, label = "p.adj.signif", 
    remove.bracket = F, tip.length = 0.01,  hide.ns = T) +
  labs(caption = title) -> psave

psave
