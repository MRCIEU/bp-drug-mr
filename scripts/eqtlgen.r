library(ieugwasr)
library(dplyr)
library(data.table)
library(biomaRt)
library(here)

# List of BP drug genes
genes <- scan(here("data","genes.txt"), "character")

# Get ensembl gene IDs
mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(mart = mart, attributes = c('hgnc_symbol','ensembl_gene_id'), filters = "hgnc_symbol", values=genes, uniqueRows = TRUE)

# Lookup in eqtlgen
ensids <- paste0("eqtl-a-", annotLookup$ensembl_gene_id)
eqtlgen <- tophits(ensids, pval=1e-6)
eqtlgen$ensembl_gene_id <- gsub("eqtl-a-", "", eqtlgen$id)
eqtlgen <- inner_join(eqtlgen, annotLookup, by="ensembl_gene_id")
eqtlgen$pop <- "EUR"
table(eqtlgen$hgnc_symbol)

# Save
save(eqtlgen, file=here("data", "eqtlgen.rdata"))


