# RICARDO GOMEZ-REYES
# 2023

# MIRGENE DB WEB SCRAPING:
library(rvest)

webpage <- read_html("https://mirgenedb.org/browse/ALL")

# Select the table using CSS selector

table_node <- html_nodes(webpage, "table")

# Extract the table content

table_content <- html_table(table_node)[[1]]

dim(table_content)

# MUST COINCIDE W/ 16667 microRNA genes ()

m <- as(table_content, "matrix")

colnames(m) <- NULL

head(m)

colNames <- m[1,]

colNames[1] <- "MirGeneDB_ID"

colNames <- gsub(" ", "_", colNames)

colnames(m) <- colNames

head(m <- m[-1,])

dim(m)

out <- m %>% as_tibble()

out <- out %>% mutate_at(vars(c("Start", "End", "UG" ,"UGUG","CNNC")), as.numeric)

wd <- "/Users/cigom/Documents/MIRNA_HALIOTIS/MIRGENEDB_20230314/"

write_tsv(out, file = paste0(wd, "/MIRGENEDB_2.1.tsv"))

# END
