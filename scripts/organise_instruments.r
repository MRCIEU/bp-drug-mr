library(here)
library(tidyverse)
library(data.table)

# Get MESA eQTLs

fn <- tibble(
    path=list.files("../../data/bp-drug-mr/extract", full.names=TRUE), 
    fn=basename(path), 
    pop=strsplit(fn, "_") %>% sapply(., function(x) x[1])
)

o <- lapply(1:nrow(fn), function(i)
{
    fread(fn$path[i]) %>%
    mutate(pop=fn$pop[i])
}) %>% bind_rows()

names(o) <- scan("../../data/bp-drug-mr/raw/AFA_cis_eqtl_summary_statistics.txt.gz", nline=1, what="character") %>% {c(., "pop")}

table(o$pop, o$gene_name)
ggplot(o, aes(x=pos_snps, y=-log10(pvalue))) +
geom_point() +
facet_grid(pop ~ gene_name, scale="free_x")


# Get ARIC pQTLs

seqid <- fread("../../data/bp-drug-mr/raw/AA/seqid.txt")
genes <- scan("genes.txt", what="character")
seqid <- subset(seqid, entrezgenesymbol %in% genes)

p <- lapply(1:nrow(seqid), function(i)
{
    fread(file.path("../../data/bp-drug-mr/raw/AA/", paste0(seqid$seqid_in_sample[i], ".PHENO1.glm.linear"))) %>% mutate(gene_name=seqid$entrezgenesymbol[i], pop="AFR")
}) %>% bind_rows()


seqid <- fread("../../data/bp-drug-mr/raw/EA/seqid.txt")
seqid <- subset(seqid, entrezgenesymbol %in% genes)

p1 <- lapply(1:nrow(seqid), function(i)
{
    fread(file.path("../../data/bp-drug-mr/raw/EA/", paste0(seqid$seqid_in_sample[i], ".PHENO1.glm.linear"))) %>% mutate(gene_name=seqid$entrezgenesymbol[i], pop="EUR", P=as.numeric(P))
}) %>% bind_rows()
str(p1)

p <- bind_rows(p, p1)
p %>% mutate(fdr=p.adjust(P)) %>%
    group_by(pop, gene_name) %>% 
    summarise(n=n(), minp=min(P), minfdr=min(fdr))

ggplot(p, aes(x=POS, y=-log10(P))) +
geom_point() +
facet_grid(pop ~ gene_name, scale="free_x")

min(o$pvalue)
o$fdr <- p.adjust(o$pvalue, "fdr")
min(o$fdr)


a <- fread("gene_info.txt")
aliases <- group_by(a, symbol) %>%
    do({
        x <- .
        alias=strsplit(x$alias, ", ") %>% unlist()
        tibble(gene_name=x$symbol, alias)

    })
aliases

seqid <- fread("../../data/bp-drug-mr/raw/EA/seqid.txt")
table(aliases$alias %in% seqid$entrezgenesymbol)
write.table(aliases$alias, file="gene_alias.txt", row=F, col=F, qu=F)


require('biomaRt')

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)

annot <- subset(annotLookup, hgnc_symbol %in% genes)
annot
ensids <- paste0("ebi-a-", annot$ensembl_gene_id)

eqtlgen <- tophits(ensids)

library(ieugwasr)
a

