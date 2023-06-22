library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(SeqVarTools)
library(qqman)


# Edit work directory

dir <- c("/user/work/ac14629/MRC_network_project/data/UKB/plink")
diro <- c("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS") 

#Read in input data
phen_bi <- fread(paste0("/user/work/ac14629/MRC_network_project/data/UKB/phenotype/mrc_phenotypes_afr.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)
cov <- fread(paste0("/user/work/ac14629/MRC_network_project/data/UKB/phenotype/mrc_covariates.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)
pc <- fread(paste0("/user/work/ac14629/MRC_network_project/data/UKB/phenotype/mrc_pca.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)

#Join phenotype, covariate and PC data files and prepare input data
cov_select <- cov %>% select(IID, age, sex) 
cov_replace <- cov_select %>% mutate(sex = ifelse(sex == 1,"M","F"))
phen_select_bi <- phen_bi %>% select("IID", "stroke") #change phenotype
dat <- merge(phen_select_bi, cov_replace, by="IID") %>% merge(., pc, by="IID") 
dat <- dat %>% na.omit()

input_dat <- data.frame(scanID = dat$ieu, pheno = dat$stroke, age = dat$age, sex = dat$sex, pc1 = dat$PC1, pc2 = dat$PC2, pc3 = dat$PC3, pc4 = dat$PC4, pc5 = dat$PC5, pc6 = dat$PC6, pc7 = dat$PC7, pc8 = dat$PC8, pc9 = dat$PC9, pc10= dat$PC10, pc11 = dat$PC11, pc12 = dat$PC12, pc13 = dat$PC13, pc14 = dat$PC14, pc15 = dat$PC15, pc16 = dat$PC16, pc17 = dat$PC17, pc18 = dat$PC18, pc19 = dat$PC19, pc20 = dat$PC20)
scanAnnot <- ScanAnnotationDataFrame(input_dat)
scanAnnot

### Make a genetic relationship matrix (GRM) from kinship estimates calculated from PC-Relate analysis
load(file = paste0(diro, "/AFR/pcrelate_kinship_AFR.RData"))
kinship_scaled <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2), scaleKin = 2)

#Fit null model 
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c("age","sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "pc11", "pc12", "pc13", "pc14", "pc15", "pc16", "pc17", "pc18", "pc19", "pc20"), cov.mat = kinship_scaled, family = "binomial")

# create a genotype object
geno <- GdsGenotypeReader(filename = paste0(dir, "/AFR/merged_chr1_22/AFR.genotype.gds"))
genoData <- GenotypeData(geno)

#Run SNP-phenotype association test
genoIterator <- GenotypeBlockIterator(genoData)
assoc <- assocTestSingle(genoIterator, null.model = nullmod, test="Score.SPA",BPPARAM = BiocParallel::SerialParam())
ea <- effectAllele(genoData, variant.id=assoc$variant.id)
assoc$chr <- as.numeric(assoc$chr)
assoc_out <- merge(assoc, ea, by="variant.id")
assoc_out_sorted <- assoc_out[order(assoc_out$chr),]

#Save output files 
write.table(assoc_out_sorted, "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/AFR/GENESIS_assoc_stroke_out.txt", row.names=F, quote=F)

#Create Manhattan and Q-Q plot 
png(file= "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/AFR/GENESIS_stroke_assoc_manhattan.png", width=8, height=7, units="in", res=300)
manhattan(assoc_out_sorted, chr="chr", bp="pos", snp="variant.id", p="SPA.pval" )
dev.off()

png(file= "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/AFR/GENESIS_stroke_assoc_qq.png", width=8, height=7, units="in", res=300)
qq(assoc_out_sorted$SPA.pval)
dev.off()
