# Create a present/absent tree of known mirna found in red abalone
# cladistics approach

DB %>% 
  mutate(Node = `Node_of_origin_(family)`) %>%
  distinct(Name, Node) %>%
  dplyr::count(Node, sort = T)

m <- DB %>%
  distinct(MirGeneDB_ID, `Node_of_origin_(family)`) %>%
  mutate(bin = 1) %>%
  pivot_wider(names_from = `Node_of_origin_(family)`, values_from = bin, values_fill = 0) %>%
  data.frame()

DB %>%
  distinct(MirGeneDB_ID, `Node_of_origin_(family)`) %>% 
  pull(`Node_of_origin_(family)`) -> labels


rownames(m) <- m$MirGeneDB_ID

m$MirGeneDB_ID <- NULL

d <- dist(m, method = "binary")

hc <- stats::hclust(d, method="ward.D2")

plot(hc)

# doi:10.1093/gbe/evy096
# We conducted phylogenetic analyses of the curated and uncurated data sets (supplementary files 1 and 2, Supplementary Material online) using a stochastic Dollo binary substitution model ... binary substitution models accommodating additional types of rate heterogeneity, as well as partial sampling schemes, may mit- igate the effect of loss, leading to greater accuracy in such data sets (Taylor et al. 2014; Thomson et al. 2014).

# remotes::install_github("liamrevell/phytools")

# https://github.com/liamrevell/Revell.phytools-v2/blob/main/Revell.phytools-v2.pdf

# 1st run Dollo assumption 

Claddis::DolloSC


library(Claddis)

# https://github.com/graemetlloyd/Claddis?tab=readme-ov-file

calculate_morphological_distances(labels)

library(phytools)

## extract discrete character (feeding mode)

m <- head(phytools::to.matrix(labels, length(labels)))

fmode<-setNames(labels, rownames(m))

## fit model

# tree: an object of class "phylo". 
t <- ape::as.phylo.hclust(hc)

er_model<-fitMk(t,fmode,model="ER", pi="fitzjohn")

## do stochastic mapping
smap<-simmap(er_model)

smm <- summary(smap)

plot(smm)

## plot posterior density on the number of changes
plot(density(smap),bty="l")

title(main="Posterior distribution of changes of each type",
  font.main=3)

ltt <- ltt(smap)
