# Library required packages
library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(qqman)

# King method is used for kinship estimation. This method is not robust against admixture.
# PC-AiR is used for population structure inference that is robust to known or cryptic relatedness
# The 1000G set unrelated part is used to ensure ancestry representative set for PCA building
# Then PC-Relate is used for accurate relatedness estimation in the presence of population structure, admixutre, and departures from Hardy-Weinberg equilibrium.
# The "PC-Air"-"PC-Relate" cycle is run twice (It can be run multiple times till it won't change any longer)
# Usually two runs are enough
##### First the King estimates are used for separating unrelated and related set (kinship (kinMat) and ancestry divergence (divMat))
##### Second run: first round PC-Relate result is used as kinship matrix, the KING matrix is still used for ancestry divergence

#Change working directory to paths to relevant input files: Phenotype, Covariates, PC scores
dir <- c("/user/work/ac14629/MRC_network_project/data/UKB/plink")
diro <- c("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS") #e.g. "your-home-directory/bp-drug-mr/outputs"

### Read in genotype data 
snpgdsBED2GDS(bed.fn = paste0(dir, "/SAS/merged_chr1_22/chr1-22_merged.bed"),
              bim.fn = paste0(dir, "/SAS/merged_chr1_22/chr1-22_merged.bim"),
              fam.fn = paste0(dir, "/SAS/merged_chr1_22/chr1-22_merged.fam"),
              out.gdsfn = paste0(dir, "/SAS/merged_chr1_22/SAS.genotype.gds"))
              
geno <- GdsGenotypeReader(filename = paste0(dir, "/SAS/merged_chr1_22/SAS.genotype.gds"))
genoData <- GenotypeData(geno)
showfile.gds(closeall=TRUE)


### PC-AiR

# set seed for LD pruning
set.seed(100)

# LD pruning
genoData <- snpgdsOpen(filename = paste0(dir, "/SAS/merged_chr1_22/SAS.genotype.gds"))
snpset <- snpgdsLDpruning(genoData, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.1), verbose=FALSE)
pruned <- unlist(snpset, use.names=FALSE)

## KING method of moment for the identity-by-descent (IBD) analysis
# Estimates the proportions of the genome at which two individuals share 0, 1, or 2 alleles IBD
# SNPRelate::KING-robust. Robust to pop. structure, not admixture
# With "KING-robust", the function would return the proportion of SNPs with zero IBS (IBS0) and kinship coefficient (kinship).
# IBS=Identity-By-State
# KINGmat$kinship is a symmetric matrix of pairwise kinship coefficients for every pair of individuals in the sample
# KING kinship estimates are negative for samples with different ancestry

KINGmat <- snpgdsIBDKING(genoData, sample.id=NULL, snp.id=NULL, autosome.only=TRUE,
    remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
    type=c("KING-robust"), family.id=NULL, num.thread=1L,
    useMatrix=FALSE, verbose=TRUE)
rownames(KINGmat$kinship) <- KINGmat$sample.id
colnames(KINGmat$kinship) <- KINGmat$sample.id

KINGmat.pair <- snpgdsIBDSelection(KINGmat)

write.table(KINGmat.pair, "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/kinship_KING_SAS.txt",quote=F, row.names=F)

pdf("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/KING_IBS0_SAS.pdf")
plot(KINGmat.pair$IBS0, KINGmat.pair$kinship, xlab="Proportion of Zero IBS",
    ylab="Estimated Kinship Coefficient (KING-robust)")
dev.off()

snpgdsClose(genoData)

# create a genotype object
geno <- GdsGenotypeReader(filename = paste0(dir, "/SAS/merged_chr1_22/SAS.genotype.gds"))
genoData <- GenotypeData(geno)

## Run PC-AiR to perform a Principal Components Analysis
# perform standard PCA on the 'unrelated subset', and predict PC values for the 'related subset'
# The reference set included in the unrelated subset, so ancestries will be better represented by the PCs
# use the KING estimates for both kinship (kinMat) and ancestry divergence (divMat), they are used for partitioning the sample into the 'unrelated' and 'related' subsets

mypcair <-pcair(genoData, kinobj= KINGmat$kinship, kin.thresh=2^(-9/2), divobj = KINGmat$kinship, div.thresh=-2^(-9/2), snp.include = pruned, verbose=T, eigen.cnt=200)

mypcair <-pcair(genoData, kinobj= KINGmat$kinship, kin.thresh=2^(-9/2), divobj = KINGmat$kinship, div.thresh=-2^(-9/2), snp.include = pruned, verbose=T, eigen.cnt=200, unrel.set = mypcair$unrels)

summary(mypcair)

varExp <- mypcair$values/sum(mypcair$values) * 100

write.table(varExp,"/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/VarianceExplainedNumber.txt", quote=F)

pdf("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/VarianceExplained_SAS.pdf")
#screeplot(mypcair)
plot(mypcair$values/sum(mypcair$values) * 100, xlab="PCs",ylab="Variance Explained (%)")
dev.off()

#plot PC-AiR PCs 
png(file= paste0(diro, "/SAS/GENESIS_pcair_pc1_pc2.png"), width=8, height=7, units="in", res=300)
plot(mypcair)
dev.off() 
png(file= paste0(diro, "/SAS/GENESIS_pcair_pc3_pc4.png"), width=8, height=7, units="in", res=300)
plot(mypcair, vx = 3, vy = 4)
dev.off() 
png(file= paste0(diro, "/SAS/GENESIS_pcair_pc5_pc6.png"), width=8, height=7, units="in", res=300)
plot(mypcair, vx = 5, vy = 6)
dev.off() 
png(file= paste0(diro, "/SAS/GENESIS_pcair_pc7_pc8.png"), width=8, height=7, units="in", res=300)
plot(mypcair, vx = 7, vy = 8)
dev.off() 
png(file= paste0(diro, "/SAS/GENESIS_pcair_pc9_pc10.png"), width=8, height=7, units="in", res=300)
plot(mypcair, vx = 9, vy = 10)
dev.off() 

pca.df <- data.frame(mypcair$vectors)
colnames(pca.df) <- c(paste("PC",1:200,sep=""))
write.table(pca.df, "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/EigenVectors_SAS.txt",sep="\t",quote=FALSE)

# Write table of related set
related_set = mypcair$rels
write.table(related_set, "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/Related_King_SAS.txt", quote=F)

#Write table of unrelated set 
unrelated_set = mypcair$unrels
write.table(unrelated_set, "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/Unrelated_King_SAS.txt", quote=F)

#save the PCA plots
pdf("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/PC_Air_SAS.pdf")
plot(mypcair)
par(mfrow=c(2,1))
plot(mypcair, vx=3, vy=4)
dev.off()

### Relatedness estimation adjusted for principal components (PC-Relate)

#Number of PCs used in the relatedness estimation 
npca <- 20

# run PC-Relate
genoData_pruned <- GenotypeBlockIterator(genoData, snpInclude=pruned)
mypcrelate <- pcrelate(genoData_pruned, mypcair$vectors[,1:npca], training.set = mypcair$unrels, BPPARAM = BiocParallel::SerialParam())

#save rdata object and write kinship coefficients to file 
save(mypcrelate, file ="/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/pcrelate_kinship_SAS.RData")
write.table(mypcrelate$kinBtwn,"/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/pcrelate_kinship_SAS.txt",row.names=F,col.names=T,sep="\t",quote=F)

# plot kinship vs k0 
pdf("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/pc-relate_IBS0_SAS_round1.pdf")
plot(mypcrelate$kinBtwn$k0, mypcrelate$kinBtwn$kin, xlab="Proportion of Zero IBS", ylab="Estimated Kinship Coefficient (KING-robust)")
dev.off()

### Make a genetic relationship matrix (GRM) from kinship estimates calculated from PC-Relate analysis 
kinship_scaled <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2), scaleKin = 2)
kinship.scaled.file <- paste("kinshipScaled_analysis_SAS_results_all_",npca,"PC.txt",sep="")
kinship.x<- kinship_scaled@x
write.table(kinship.x,file="/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/SAS/kinship.scaled.SAS.file",row.names=F,col.names=T,sep="\t",quote=F)
save.image("relatedness.RData")
