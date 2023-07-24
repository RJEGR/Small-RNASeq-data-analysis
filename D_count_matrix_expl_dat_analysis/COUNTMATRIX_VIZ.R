
# THIS IS THE UPGRADE VERSION PREV. TO DE ANALYSIS


rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "~/Documents/MIRNA_HALIOTIS/SHORTSTACKS/ShortStack_20230315_out/"

# 
pHpalette <- c(`Experimental`="#ad1f1f", `Control`= "#4575b4")
#

out <- read_rds(paste0(path, "Counts_and_mtd.rds"))

dim(COUNTS <- out[[1]])

dim(MTD <- out[[2]])

# PREVALENCE BY REP ----

nrow(raw_count <- COUNTS[colNames]) # 66,226

head(raw_count <- as.data.frame(raw_count, row.names = COUNTS$Name))

rownames(raw_count) <- COUNTS$Name

apply(raw_count, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(raw_count, 1, function(x) sum(x > 0))

# mean_se = apply(raw_count, 1, function(x) mean_se(x)) %>% do.call(rbind, .)

data.frame(Prevalence = prevelancedf, 
  TotalAbundance = rowSums(raw_count)) %>% # mean_se
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf <- prevelancedf %>% left_join(RESULTS %>% select(Name, MIRNA, KnownRNAs, Strand, DicerCall))

prevelancedf %>% 
  count(Prevalence) %>% 
  ggplot(aes(Prevalence, n)) + geom_col() +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Number of sRNAs", labels = scales::comma) +
  scale_x_continuous(breaks = 1:12) -> ps

ggsave(ps, filename = 'prevalence_hist.png', 
  path = path, width = 3, height = 2, device = png)

prevelancedf %>% 
  arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence))) %>%
  mutate(TotalAbundance = log2(TotalAbundance+1)) %>% # edgeR::cpm(count) 
  ggplot(aes(TotalAbundance)) + geom_histogram() + 
  facet_wrap(~ Prevalence, scales = 'free_y') -> p1

dat_text <- prevelancedf %>% group_by(Prevalence) %>% tally() %>% 
  mutate(cumsum = cumsum(n)) %>% arrange(Prevalence) %>%
  mutate(Prevalence = paste0(Prevalence, ' Samples')) %>%
  mutate(Prevalence = factor(Prevalence, levels = unique(Prevalence)))

p1 <- p1 + geom_text(
  data    = dat_text, family = "GillSans",
  mapping = aes(x = -Inf, y = -Inf, label = paste0(n, " sRNAs")),
  hjust   = -1, vjust   = -2.5, size = 2.5) + 
  theme_classic(base_size = 7, base_family = "GillSans") +
  labs(x = expression(~Log[2]~('TotalAbundance'~+1)), y = "") +
  scale_y_continuous("Number of sRNAs", labels = scales::comma)


ggsave(p1, filename = 'prevalence_hist_facet.png', 
  path = path, width = 7, height = 4, device = png)


prevelancedf %>% 
  mutate(MIRNA = factor(MIRNA, levels = c("Y","N") )) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(TotalAbundance, Prevalence/12)) +
  # geom_point(aes(color = MIRNA), alpha = 0.2) +
  stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE) +
  scale_x_log10("Total Abundance (log10 scale)", labels = scales::comma) +  
  scale_y_continuous("Prevalence [Frac. Samples]", 
    labels = scales::percent_format(scale = 100)) +
  theme_classic(base_family = "GillSans") +
  theme(legend.position="top") -> ps
# facet_grid(~ MIRNA) 
# scale_color_manual(values = c("red", "grey"))


ggsave(ps, filename = 'prevalence.png', 
  path = path, width = 4, height = 4, device = png)


prevelancedf %>% 
  mutate(Prevalence = factor(Prevalence)) %>%
  # filter(MIRNA == "Y") %>%
  ggplot(aes(log10(TotalAbundance), color  = Prevalence, fill = Prevalence)) +
  # scale_x_continuous() +
  stat_ecdf() 
# geom_density(alpha = 0.3)



