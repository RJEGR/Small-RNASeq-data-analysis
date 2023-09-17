library(AnnotationHub)
package = "AnnotationHub"

oldcache = path.expand(rappdirs::user_cache_dir(appname=package))
setAnnotationHubOption("CACHE", oldcache)
ah = AnnotationHub(localHub=TRUE)
## removes old location and all resources
removeCache(ah, ask=FALSE)

## create the new default caching location
newcache = tools::R_user_dir(package, which="cache")
setAnnotationHubOption("CACHE", newcache)
ah = AnnotationHub()

unique(ah$species)[grepl("Haliotis ruf", unique(ah$species))]

# query(ah, "OrgDb")

library(AnnotationHub)

q <- query(ah, "Haliotis")
id <- q$ah_id[1]

Haliotis <- ah[[id]]

keytypes(Haliotis)

columns(Haliotis)


# OrgDb <- query(ah, "OrgDb")

# unique(OrgDb$species)[grepl("Haliotis", unique(OrgDb$species))]

# OrgDb <- query(ah, "OrgDb")[which(grepl("Haliotis rufescens", unique(OrgDb$species)))]

GOSemSim::godata(Haliotis, ont="BP")

GOSemSim::load_OrgDb(Haliotis)

kk <- keys(Haliotis, keytype = "ENTREZID")

goAnno <- suppressMessages(select(Haliotis, keys = kk, keytype = "ENTREZID", 
  columns = columns(Haliotis)))

# c("GID", "ONTOLOGY")

dim(goAnno)

# NOT GO IN ORGDB BYT ENTREZID, SO , LETS TRY TO USE GETGO W/ BIOMART

ensembl_metazoa <- useEnsembl(biomart = "metazoa_mart",
  host = "https://metazoa.ensembl.org")

datasets <- listDatasets(ensembl_metazoa)

# head(datasets)

datasets <- datasets %>% filter(str_detect(description, 'Haliotis rufescens'))

# searchDatasets(mart = ensembl_metazoa, pattern = "cgigas")

ensembl_db <- useEnsembl(
  biomart = "metazoa_mart",
  host = "https://metazoa.ensembl.org",
  dataset = datasets$dataset)


mart <- biomaRt::useEnsembl("ensembl")

GETGO <- getEN(go, de_df, pattern, GOtype, 
  ensembl_mart)

GETGO <- getBM(attributes = c('entrezgene_id', 
  'go_id',
  'name_1006'),
  filters = 'go', 
  values = GO_universe, 
  mart = mart)

# ensembl_mart <- biomaRt::useEnsembl(biomart = "metazoa_mart",
#   host = "metazoa.ensembl.org",
#   dataset = 'cgigas_eg_gene')