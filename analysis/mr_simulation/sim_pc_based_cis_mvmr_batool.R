###Simulate regional summary data with LD matrix###

#First get LD matrix for a region for each ancestry

library(simulateGP)
library(dplyr)
library(MendelianRandomization)
library(stats)

load(url("https://github.com/explodecomputer/simulateGP/raw/master/data/ldetect.rdata"))
head(ldetect)
a <- subset(ldetect, pop == "EUR") %>% mutate(len=stop-start) %>% arrange(len) %>% slice(100)

eurbfile <- "/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EUR"
afrbfile <- "/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/AFR"
easbfile <- "/Users/gh13047/repo/opengwas-api-internal/opengwas-api/app/ld_files/EAS"

eurld <- get_ld(a$chr, a$start, a$stop, eurbfile)
afrld <- get_ld(a$chr, a$start, a$stop, afrbfile)
easld <- get_ld(a$chr, a$start, a$stop, easbfile)

#Then simulate summary data
#
#- 10 causal variants for trait 1
#- Those variants have an effect on trait 2 + 3 other variants
#- Those variants have an effect on trait 3 + 3 other variants
#- Those variants have an effect on outcome + 3 other variants

rsids <- Reduce(intersect, list(eurld$map$snp, afrld$map$snp, easld$map$snp))
causalx1 <- sample(rsids, 10)
causalx2 <- sample(rsids, 3)
causalx3 <- sample(rsids, 3)
causaly <- sample(rsids, 3)

paramseurx1 <- eurld$map %>%
    generate_gwas_params(h2=0, S=0, Pi=0) %>%
    mutate(beta = case_when(snp %in% causalx1 ~ 0.05, TRUE ~ 0))

paramseurx2 <- paramseurx1 %>% mutate(beta = beta * -0.7) %>% 
    mutate(beta = case_when(snp %in% causalx2 ~ beta + 0.05, TRUE ~ beta))

paramseurx3 <- paramseurx1 %>% mutate(beta = beta * 0.5) %>% 
    mutate(beta = case_when(snp %in% causalx3 ~ beta + 0.05, TRUE ~ beta))

paramseury <- paramseurx1 %>% mutate(beta = beta * 0.5) %>% 
    mutate(beta = case_when(snp %in% causaly ~ beta + 0.05, TRUE ~ beta))

sseurx1 <- paramseurx1 %>% generate_gwas_ss(nid=100000, ldobj=eurld)
sseurx2 <- paramseurx2 %>% generate_gwas_ss(nid=100000, ldobj=eurld)
sseurx3 <- paramseurx3 %>% generate_gwas_ss(nid=100000, ldobj=eurld)
sseury <- paramseury %>% generate_gwas_ss(nid=100000, ldobj=eurld)


## Perform analysis using PCs using Batool et al. 2022 approach

# Choose a region with <=4000 SNPs otherwise it gets very slow
# Make sure that there is complete overlap of SNPs across the traits (e.g. remove SNPs that aren't present in all)
# Do this procedure for each ancestry separately

# Code below
# sseurx1 = summary stats for trait 1 (e.g. LDL)
# sseurx2 = summary stats for trait 2 (e.g. HDL)
# sseurx3 = summary stats for trait 3 (e.g. Trigs)
# sseury = summary stats for outcome (e.g. CHD)

# bhat = effect estimate
# se = standard error

# The LD matrix should be harmonised with the data, e.g. same set of SNPs, in the same order, and with the same alleles
# e.g. to make the alleles the same create a vector 'flip' that is 1 if the LD + allele is the same as the beta effect allele, and 0 if not, then
# m <- flip %*% t(flip)
# ld <- ld * m

# MV-PCA method (Batool et al. 2022)
Psi = ((abs(sseurx1$bhat)+abs(sseurx2$bhat)+abs(sseurx3$bhat))/sseury$se)%o%
((abs(sseurx1$bhat)+abs(sseurx2$bhat)+abs(sseurx3$bhat))/sseury$se)*eurld$ld

K = which(cumsum(prcomp(Psi, scale=FALSE)$sdev^2/
sum((prcomp(Psi, scale=FALSE)$sdev^2)))>0.99)[1]
K

eurx1comp = as.numeric(sseurx1$bhat%*%prcomp(Psi, scale=FALSE)$rotation[,1:K])
eurx2comp = as.numeric(sseurx2$bhat%*%prcomp(Psi, scale=FALSE)$rotation[,1:K])
eurx3comp = as.numeric(sseurx3$bhat%*%prcomp(Psi, scale=FALSE)$rotation[,1:K])
eurycomp = as.numeric(sseury$bhat%*%prcomp(Psi, scale=FALSE)$rotation[,1:K])
eurycomp_se = as.numeric(sseury$se%*%prcomp(Psi, scale=FALSE)$rotation[,1:K])

# To do weights
Omega <- sseury$se %o% sseury$se * eurld$ld
pcOmega <- t(prcomp(Psi, scale=FALSE)$rotation[,1:K])%*%Omega%*%
prcomp(Psi, scale=FALSE)$rotation[,1:K]

# Univariable MR (weighted)
summary(lm(eurycomp ~ 0 + eurx1comp, weight = 1/diag(pcOmega)^2))

# MV-PCA (weighted)
mvpca = mr_mvivw(mr_mvinput(cbind(eurx1comp, eurx2comp, eurx3comp),
cbind(rep(1, length(eurx1comp)), rep(1, length(eurx1comp)), rep(1, length(eurx1comp))),
eurycomp, rep(1, length(eurx1comp)), corr=pcOmega), model="fixed")
mvpca$Estimate; mvpca$StdError

mvpca_est = solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%
cbind(eurx1comp, eurx2comp, eurx3comp))%*%
rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%eurycomp
mvpca_se = sqrt(diag(solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%
cbind(eurx1comp, eurx2comp, eurx3comp))))

z_scores <- mvpca_est / mvpca_se
p_values <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
p_values
