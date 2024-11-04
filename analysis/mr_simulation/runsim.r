library(simulateGP)
library(dplyr)
library(data.table)
library(stats)
library(MVMR)
library(furrr)
library(here)
library(ggplot2)


sim_dat <- function(ld, ncausal, nid, bgx=0.05) {

    rsids <- ld$map$snp
    causalx1 <- sample(rsids, ncausal)
    causalx2 <- sample(rsids, ncausal)
    causalx3 <- sample(rsids, ncausal)
    causaly <- sample(rsids, ncausal)

    paramseurx1 <- ld$map %>%
        generate_gwas_params(h2=0, S=0, Pi=0) %>%
        mutate(beta = case_when(snp %in% causalx1 ~ bgx, TRUE ~ 0))

    paramseurx2 <- paramseurx1 %>% mutate(beta = beta * -0.7) %>% 
        mutate(beta = case_when(snp %in% causalx2 ~ beta + bgx, TRUE ~ beta))

    paramseurx3 <- paramseurx1 %>% mutate(beta = beta * 0.5) %>% 
        mutate(beta = case_when(snp %in% causalx3 ~ beta + bgx, TRUE ~ beta))

    paramseury <- paramseurx1 %>% mutate(beta = beta * 0.5) %>% 
        mutate(beta = case_when(snp %in% causaly ~ beta + bgx, TRUE ~ beta))

    sseurx1 <- paramseurx1 %>% generate_gwas_ss(nid=nid, ldobj=ld)
    sseurx2 <- paramseurx2 %>% generate_gwas_ss(nid=nid, ldobj=ld)
    sseurx3 <- paramseurx3 %>% generate_gwas_ss(nid=nid, ldobj=ld)
    sseury <- paramseury %>% generate_gwas_ss(nid=nid, ldobj=ld)

    return(list(sseurx1=sseurx1, sseurx2=sseurx2, sseurx3=sseurx3, sseury=sseury))
}


sim_est <- function(sseurx1, sseurx2, sseurx3, sseury, ld) {

    Psi = ((abs(sseurx1$bhat) + abs(sseurx2$bhat) + abs(sseurx3$bhat)) / sseury$se) %o% 
        ((abs(sseurx1$bhat) + abs(sseurx2$bhat) + abs(sseurx3$bhat)) / sseury$se) * 
        ld$ld
    pcs <- eigen(Psi)

    K <- which(cumsum(pcs$values^2) / sum(pcs$values^2) > 0.99)[1]

    eurx1comp <- sseurx1$bhat %*% pcs$vectors[,1:K] %>% drop()
    eurx2comp <- sseurx2$bhat %*% pcs$vectors[,1:K] %>% drop()
    eurx3comp <- sseurx3$bhat %*% pcs$vectors[,1:K] %>% drop()
    eurycomp <- sseury$bhat %*% pcs$vectors[,1:K] %>% drop()
    eurycomp_se <- sseury$se %*% pcs$vectors[,1:K] %>% drop()

    dat <- tibble(eurx1comp, eurx2comp, eurx3comp, eurycomp, eurycomp_se, eurx1comp_se = 1, eurx2comp_se = 1, eurx3comp_se = 1)

    d <- format_mvmr(
        cbind(dat$eurx1comp, dat$eurx2comp, dat$eurx3comp),
        dat$eurycomp,
        cbind(dat$eurx1comp_se, dat$eurx2comp_se, dat$eurx3comp_se),
        dat$eurycomp_se,
        1:K
    )
    f <- suppressWarnings(strength_mvmr(d)) %>% unlist()



    Omega <- sseury$se %o% sseury$se * ld$ld
    pcOmega <- t(pcs$vectors[,1:K]) %*% Omega %*% pcs$vectors[,1:K]


    mvpca_est = solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%
    cbind(eurx1comp, eurx2comp, eurx3comp))%*%
    rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%eurycomp
    mvpca_se = sqrt(diag(solve(rbind(eurx1comp, eurx2comp, eurx3comp)%*%solve(pcOmega)%*%
    cbind(eurx1comp, eurx2comp, eurx3comp))))

    z_scores <- mvpca_est / mvpca_se
    p_values <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
    p_values

    res <- tibble(x=paste0("x", 1:3), b=mvpca_est %>% drop(), se=mvpca_se %>% drop(), z=z_scores %>% drop(), p=p_values %>% drop(), K=K, Fstat=f)
    return(res)
}


load(url("https://github.com/explodecomputer/simulateGP/raw/master/data/ldetect.rdata"))
head(ldetect)
a <- subset(ldetect, pop == "EUR") %>% mutate(len=stop-start) %>% arrange(len) %>% slice(100)

eurbfile <- "/local-scratch/data/1000g/EUR"
eurld <- get_ld(a$chr, a$start, a$stop, eurbfile)

str(eurld)


ld1 <- eurld
ld1$ld <- ld1$ld[1:1000, 1:1000]
ld1$map <- ld1$map[1:1000,]
str(ld1)

ld2 <- eurld
ld2$ld <- ld2$ld[1001:2000, 1001:2000]
ld2$map <- ld2$map[1001:2000,]
str(ld2)

ld3 <- eurld
ld3$ld <- ld3$ld[2001:3000, 2001:3000]
ld3$map <- ld3$map[2001:3000,]
str(ld3)

ld4 <- eurld
ld4$ld <- ld4$ld[3001:4000, 3001:4000]
ld4$map <- ld4$map[3001:4000,]
str(ld4)


d <- sim_dat(ld1, 3, 100000)

sim_est(d$sseurx1, d$sseurx2, d$sseurx3, d$sseury, ld1)


sim <- function(ld_name, nid, ncausal, bgx=0.05, sim=1) {
    args <- as.list(environment()) %>% as_tibble()
    ld <- get(ld_name)
    d <- sim_dat(ld, ncausal, nid)
    e <- sim_est(d$sseurx1, d$sseurx2, d$sseurx3, d$sseury, ld)
    e <- bind_cols(e, args)
    return(e)
}

param <- expand.grid(
    ld_name=c("ld1", "ld2", "ld3", "ld4"),
    nid=c(10000, 100000, 1000000),
    ncausal=c(3, 5, 10),
    bgx=c(0.05, 0.1, 0.2),
    sim = 1:20,
    stringsAsFactors=FALSE
)

dim(param)


a <- sim("ld1", 100000, 10, 0.2)
a

plan(multicore, workers = 100)
op <- furrr_options(seed=TRUE)
res <- future_pmap(param, sim, .progress=TRUE, .options=op) %>% bind_rows()

saveRDS(res, file=here("analysis/mr_simulation/results/", "mrsims.rds"))



res <- readRDS(here("analysis/mr_simulation/results/", "mrsims.rds"))


p <- ggplot(res, aes(x=as.factor(x), y=b)) +
geom_boxplot(outliers=FALSE)
ggsave(p, file=here("analysis/mr_simulation/results/", "temp1.png"))

p <- ggplot(res, aes(x=as.factor(x), y=z)) +
geom_boxplot(outliers=FALSE)
ggsave(p, file=here("analysis/mr_simulation/results/", "temp2.png"))

p <- ggplot(res, aes(x=as.factor(x), y=b)) +
geom_boxplot(outliers=FALSE, aes(fill=as.factor(nid))) +
facet_wrap(bgx ~ ncausal) +
geom_hline(yintercept=0, linetype="dashed", color="red") +
geom_hline(yintercept=0.5, linetype="dashed", color="black")
ggsave(p, file=here("analysis/mr_simulation/results/", "temp3.png"))


p <- ggplot(res, aes(x=as.factor(x), y=Fstat)) +
geom_boxplot(outliers=FALSE, aes(fill=as.factor(nid))) +
facet_wrap(bgx ~ ncausal) +
geom_hline(yintercept=0, linetype="dashed", color="red") +
geom_hline(yintercept=0.5, linetype="dashed", color="black")
ggsave(p, file=here("analysis/mr_simulation/results/", "temp4.png"))


res$btrue <- 0
res$btrue[res$x == "x1"] <- 0.5

p <- ggplot(res, aes(x=as.factor(x), y=(b-btrue)^2)) +
geom_boxplot(outliers=FALSE, aes(fill=as.factor(nid))) +
facet_wrap(bgx ~ ncausal)
ggsave(p, file=here("analysis/mr_simulation/results/", "temp5.png"))

res$mse <- (res$b - res$btrue)^2

summary(lm(mse ~ nid + ncausal + bgx + x, data=res))
summary(lm(abs(z) ~ nid + ncausal + bgx, data=res %>% filter(x=="x1")))

