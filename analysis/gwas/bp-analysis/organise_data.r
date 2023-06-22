# Outline
# This script extracts the Gill defined SBP drug target regions from each population and harmonises them

# 1. Use multi-ancestry meta-analysis to identify ancestry-agnostic top-hits (need to add chr:pos from dbSNP VCF)
# 2. Use Gill regions to extract best SNP from each potential region
# 3. Use UKBB cross-ancestry panel to clump
# 4. Standardise and harmonise across ancestries using chr:pos_a1_a2 (a1,a2 alphabetical)
# 5. Save extracted data


library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(ieugwasr)
library(parallel)

set_bcftools()
bcftools <- options()[["tools_bcftools"]]

# Read in dataset
dat <- fread("/projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/GWAS/Meta_analyses/METAL_METAANALYSIS_SBP_UKB_UGR.TBL")  

# Need to get chr:pos
# Use dbsnp vcf file to lookup all rsids and then extract chr, pos, ID via bcftools
rsid <- dat$MarkerName
vcffile <- "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/dbsnp/dbsnp.v153.b37.vcf.gz"
tmp <- tempfile()
write.table(unique(rsid), file=paste0(tmp, ".snplist"), row.names = FALSE, col.names = FALSE, quote = FALSE)
cmd <- glue("{bcftools} query -f '%CHROM %POS %ID\\n' -i'ID=@{tmp}.snplist' {vcffile} > {tmp}.txt")
system(cmd)

map <- fread(glue("{tmp}.txt"))
dat <- left_join(dat, map, by=c("MarkerName" = "V3"))


# Sort out non rsid IDs
chrpos <- grep("rs", rsid, value=TRUE, invert=TRUE)
a <- do.call(rbind, strsplit(chrpos, split=":"))
b <- sapply(strsplit(a[,2], split="_"), \(x) x[1])
a[,2] <- b
chrpos <- tibble(MarkerName=chrpos, V1=a[,1], V2=a[,2])
index <- match(chrpos$MarkerName, dat$MarkerName)
stopifnot(all(dat$MarkerName[index] == chrpos$MarkerName))
dat$V1[index] <- chrpos$V1
dat$V2[index] <- chrpos$V2
subset(dat, MarkerName == chrpos$MarkerName[100])
dat$V2 <- as.numeric(dat$V2)

table(is.na(dat$V2))
table(is.na(dat$V1))
subset(dat, is.na(V1))

# Extract Gill chr:pos

gill <- fread("data/gill_chrpos.csv") %>%
    as_tibble() %>%
tidyr::separate(`Position (hg19)`, sep="-", into=c("start", "end"))
gill

dat2 <- subset(dat, `P-value` < 1e-4 & !grepl("\\?", Direction))
dim(dat2)
extr <- mclapply(1:nrow(gill), \(i)
{
    a <- subset(dat2, V1 == gill$Chromosome[i] & V2 >= gill$start[i] & V2 <= gill$end[i])
    a$`Target gene` <- gill$`Target gene`[i]
    a$`Function` <- gill$`Function`[i]
    a$`Drug` <- gill$`Drug`[i]
    a
}, mc.cores=10) %>% bind_rows()

# Clump using European reference panel
extr$rsid <- extr$MarkerName
extr$pval <- extr$`P-value`
bfile <- "/projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/for_gib/LD_ref/chr1-22_merged"
plink_bin <- "plink"
extr_clumped <- ld_clump(extr, plink_bin=plink_bin, bfile=bfile) %>%
    filter(!duplicated(rsid))

dim(extr_clumped)
table(extr_clumped$`Drug`)
table(gill$`Drug`)

# Write
write.csv(extr_clumped, file="data/meta_gill_extract.csv")

# Need to extract relevant SNPs from each of the ancestry-specific GWASs

clumped <- fread("/projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/GWAS/Meta_analyses/METAANALYSIS_SBP_UKB_UGR_clumped")

sel_rsid <- unique(c(clumped$rsid, extr_clumped$rsid))
tmp <- tempfile()
write.table(sel_rsid, file=tmp, row=F, col=F, qu=F)

pops <- c("afr", "eas", "eur", "sas", "ugr")

extract <- mclapply(pops, \(i) {
    cmd <- glue("grep -wf {tmp} /projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/for_gib/{i}_sbp_assoc_imputed.txt > {tmp}_{i}")
    system(cmd)
    fread(glue("{tmp}_{i}"))
}, mc.cores=6)

ex <- list()

a <- scan("/projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/for_gib/afr_sbp_assoc_imputed.txt", what="character", nlines=1)
names(extract[[1]]) <- a
ex$afr <- extract[[1]] %>% 
    as_tibble() %>%
    dplyr::select(rsid=SNP, chr, pos, eaf=freq, ea=effect.allele, oa=other.allele, pval=Score.pval, beta=Est, se=Est.SE, n=n.obs) %>% 
    mutate(pop="afr")

names(extract[[2]]) <- a
ex$eas <- extract[[2]] %>% 
    as_tibble() %>%
    dplyr::select(rsid=SNP, chr, pos, eaf=freq, ea=effect.allele, oa=other.allele, pval=Score.pval, beta=Est, se=Est.SE, n=n.obs) %>% 
    mutate(pop="eas")

names(extract[[4]]) <- a
ex$sas <- extract[[4]] %>% 
    as_tibble() %>%
    dplyr::select(rsid=SNP, chr, pos, eaf=freq, ea=effect.allele, oa=other.allele, pval=Score.pval, beta=Est, se=Est.SE, n=n.obs) %>% 
    mutate(pop="sas")

names(extract[[5]]) <- a
ex$ugr <- extract[[5]] %>% 
    as_tibble() %>%
    dplyr::select(rsid=SNP, chr, pos, eaf=freq, ea=effect.allele, oa=other.allele, pval=Score.pval, beta=Est, se=Est.SE, n=n.obs) %>% 
    mutate(pop="ugr")

a <- scan("/projects/MRC-IEU/research/projects/ieu2/p4/080/working/results/for_gib/eur_sbp_assoc_imputed.txt", what="character", nlines=1)
names(extract[[3]]) <- a
ex$eur <- extract[[3]] %>%
    as_tibble() %>%
    dplyr::select(rsid=SNP, chr=CHR, pos=BP, eaf=A1FREQ, ea=ALLELE1, oa=ALLELE0, pval=P_BOLT_LMM, beta=BETA, se=SE) %>% 
    mutate(pop="eur")

standardise <- function(d)
{
    toflip <- d$ea > d$oa
    d$eaf[toflip] <- 1 - d$eaf[toflip]
    d$beta[toflip] <- d$beta[toflip] * -1
    temp <- d$oa[toflip]
    d$oa[toflip] <- d$ea[toflip]
    d$ea[toflip] <- temp
    d$snpid <- paste0(d$chr, ":", d$pos, "_", d$ea, "_", d$oa)
    d %>% arrange(chr, pos)
}

ex <- lapply(ex, standardise) %>%
    bind_rows()

ex <- left_join(ex, extr_clumped %>% dplyr::select(rsid, Drug, `Target gene`, Function), by="rsid")

saveRDS(ex, file="data/pop_extract_clumped.rds")

