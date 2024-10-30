# We replaced correlative models of targeting efficacy with a principled, biochem- ical model that explains and predicts about half of the variability attributable to the direct effects of miRNAs on their targets (McGeary et al 2019). 

# ScanMir https://doi.org/10.1093/bioinformatics/btac110

# LOAD SET OF 000 genes
# Use ID to filter only 3 UTR 
# Load 118 major RNAs
# Run model to calculate targeting efficacy
# Para construir tu modelo debes remplazar el objeto dummyKdData() por una tabla de valores kd reportados para los microRNAs que identifique en el abulon

rm(list = ls());

if(!is.null(dev.list())) dev.off()


dir <- "~/Documents/MIRNA_HALIOTIS/MIRNA_PUB_2024/"

print(DB <- read_tsv(file.path(dir, "LONGER_RELATIONAL_DB.tsv")))

# Get UTR

q.genes <- DB %>% distinct(gene_id) %>% pull(gene_id)

f <- list.files(path = "/Users/cigom/Documents/MIRNA_HALIOTIS/ENSEMBLE/", pattern = "three_prime_utr.fa", full.names = T)

seqs <- Biostrings::readDNAStringSet(f)

str(keep <- names(seqs))

str(keep <- sapply(strsplit(keep, " "), `[`, 2))

sum(keep <- keep %in% q.genes)

seqs <- seqs[keep]

# Try seed 
x <- DB %>% distinct(MajorRNA) %>% mutate(seed = substr(MajorRNA, 2,13)) %>% pull(seed, name = MajorRNA)

x <- gsub("U", "T", x)

# To predict the dissociation constant (and binding type, if any) of a given 12-mer sequence, you can use the assignKdType function:
# In this case findSeedMatches also returns the predicted affinity value for each match. The log_kd column contains log(Kd) values multiplied by 1000, where Kd is the predicted dissociation constant of miRNA:mRNA binding for the putative binding site.

data("SampleKdModel")

assignKdType(x, SampleKdModel) 

matches1 <- findSeedMatches(seqs, SampleKdModel, verbose = T)

matches1 %>% as_tibble() %>% ggplot(aes(log_kd)) + geom_histogram()

matches1 %>% as_tibble() %>% count(type)


matches1 %>% as_tibble() %>% 
  ggplot(aes(log_kd, fill = type), alpha = 0.5) + 
  geom_histogram() +
  facet_grid(type~.)

# The log_kd column contains log(Kd) values multiplied by 1000 and stored as an integer (which is more economical when dealing with millions of sites). In the example above, -5129 means -5.129, or a dissociation constant of 0.0059225. The smaller the values, the stronger the relative affinity.

# Using own model

kd <- dummyKdData()

mod3 <- getKdModel(kd=kd, mirseq="TTAATGCTAATCGTGATAGGGGTT", name = "my-miRNA")

findSeedMatches(seqs, mod3, verbose = T)

x <- DB %>% distinct(MajorRNA) %>% pull(MajorRNA)

x <- gsub("U", "T", x)

getKdModel(kd=kd, mirseq="TTATTGCTTGA", name = "my-miRNA")


# Find seed matches

y <- DB %>% distinct(MajorRNA) %>% pull(MajorRNA)

matches <- findSeedMatches(seqs, y, verbose = T)

viewTargetAlignment(matches[1], y[1], seqs)

viewTargetAlignment(matches[1], y[1], seqs)

exit



assignKdType(c("CTAGCATTAAGT","ACGTACGTACGT"), SampleKdModel)  

