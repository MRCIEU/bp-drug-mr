library(here)
library(dplyr)
library(ggplot2)
library(data.table)

# Get MESA eQTLs

fn <- tibble(
    path=list.files(here("data/extract"), full.names=TRUE), 
    fn=basename(path), 
    pop=strsplit(fn, "_") %>% sapply(., function(x) x[1])
)

o <- lapply(1:nrow(fn), function(i)
{
    fread(fn$path[i]) %>%
    mutate(pop=fn$pop[i])
}) %>% bind_rows()

names(o) <- scan(here("data/raw/AFA_cis_eqtl_summary_statistics.txt.gz"), nline=1, what="character") %>% {c(., "pop")}

table(o$pop, o$gene_name)

ggplot(o, aes(x=pos_snps, y=-log10(pvalue))) +
geom_point() +
facet_grid(pop ~ gene_name, scale="free_x")

o$fdr <- p.adjust(o$pvalue, "fdr")
min(o$fdr)

# Get ARIC pQTLs

genes <- scan(here("data", "genes.txt"), what="character")
seqid <- fread(here("data/raw/AA/seqid.txt"))
seqid <- subset(seqid, entrezgenesymbol %in% genes)

p1 <- lapply(1:nrow(seqid), function(i)
{
    fread(here("data/raw/AA/", paste0(seqid$seqid_in_sample[i], ".PHENO1.glm.linear"))) %>% mutate(gene_name=seqid$entrezgenesymbol[i], pop="AFR")
}) %>% bind_rows()


seqid <- fread(here("data/raw/EA/seqid.txt"))
seqid <- subset(seqid, entrezgenesymbol %in% genes)

p2 <- lapply(1:nrow(seqid), function(i)
{
    fread(here("data/raw/EA/", paste0(seqid$seqid_in_sample[i], ".PHENO1.glm.linear"))) %>% mutate(gene_name=seqid$entrezgenesymbol[i], pop="EUR", P=as.numeric(P))
}) %>% bind_rows()

p <- bind_rows(p1, p2)
p %>% mutate(fdr=p.adjust(P)) %>%
    group_by(pop, gene_name) %>% 
    summarise(n=n(), minp=min(P), minfdr=min(fdr))

ggplot(p, aes(x=POS, y=-log10(P))) +
geom_point() +
facet_grid(pop ~ gene_name, scale="free_x")

min(o$pvalue)
o$fdr <- p.adjust(o$pvalue, "fdr")
min(o$fdr)

aric_pqtl <- p
mesa_eqtl <- o
save(aric_pqtl, file=here("data", "aric_pqtl.rdata"))
save(mesa_eqtl, file=here("data", "mesa_eqtl.rdata"))


