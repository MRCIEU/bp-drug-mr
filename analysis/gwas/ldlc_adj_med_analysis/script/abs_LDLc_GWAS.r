# GWAS of LDL cholesterol adjusted for different absolute medication corrections

# Read in relevant packages
library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(SeqVarTools)
library(qqman)

# Edit work directory
dir <- c("/user/work/ac14629/MRC_network_project/data/UKB/phenotype")
dirg <- c("/user/work/ac14629/MRC_network_project/data/UKB/plink")
dirk <- c("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS")
diro <- c("/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/ldl_analysis/")

# Read in input data
phen_con <- fread(paste0(dir,"/mrc_ldlc_adj_abs_med_afr.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)
cov <- fread(paste0(dir,"/mrc_covariates.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)
pc <- fread(paste0(dir,"/mrc_pca.txt"), header=T, fill = T, na.strings = c("", NA), data.table = F, stringsAsFactors = FALSE)

# Make a genetic relationship matrix (GRM) from kinship estimates calculated from PC-Relate analysis
load(file = paste0(dirk, "/AFR/pcrelate_kinship_AFR.RData"))
kinship_scaled <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2), scaleKin = 2)

# Reformat covariate data
cov_select <- cov %>% select(IID, age, sex)
cov_replace <- cov_select %>% mutate(sex = ifelse(sex == 1, "M", "F"))

#Close gds file
showfile.gds(closeall=TRUE)

# Define a vector of phenotype names
phenotype_names <- c("ldlc_0.3", "ldlc_0.5", "ldlc_0.7", "ldlc_0.9", "ldlc_1.0", "ldlc_1.2", "ldlc_1.4", "ldlc_1.6", "ldlc_1.8", "ldlc_2.0")

# Loop over each phenotype
for (phenotype in phenotype_names) {
  phen_select_con <- data.frame(IID = phen_con$IID, pheno = phen_con[[phenotype]])

  dat <- merge(phen_select_con, cov_replace, by = "IID", all.x = TRUE) %>% merge(., pc, by = "IID", all.x = TRUE)
  dat <- dat %>% na.omit()

  input_dat <- data.frame(
    scanID = dat$ieu,
    pheno = dat$pheno,
    age = dat$age,
    sex = dat$sex,
    pc1 = dat$PC1,
    pc2 = dat$PC2,
    pc3 = dat$PC3,
    pc4 = dat$PC4,
    pc5 = dat$PC5,
    pc6 = dat$PC6,
    pc7 = dat$PC7,
    pc8 = dat$PC8,
    pc9 = dat$PC9,
    pc10 = dat$PC10,
    pc11 = dat$PC11,
    pc12 = dat$PC12,
    pc13 = dat$PC13,
    pc14 = dat$PC14,
    pc15 = dat$PC15,
    pc16 = dat$PC16,
    pc17 = dat$PC17,
    pc18 = dat$PC18,
    pc19 = dat$PC19,
    pc20 = dat$PC20
  )

  scanAnnot <- ScanAnnotationDataFrame(input_dat)

  # Fit null model
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c("age", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10", "pc11", "pc12", "pc13", "pc14", "pc15", "pc16", "pc17", "pc18", "pc19", "pc20"), cov.mat = kinship_scaled, family = "gaussian")

  # create a genotype object
  geno <- GdsGenotypeReader(filename = paste0(dirg, "/AFR/merged_chr1_22/AFR.genotype.gds"))
  genoData <- GenotypeData(geno)

  # Run SNP-phenotype association test
  genoIterator <- GenotypeBlockIterator(genoData)
  assoc <- assocTestSingle(genoIterator, null.model = nullmod, test = "Score", BPPARAM = BiocParallel::SerialParam())
  ea <- effectAllele(genoData, variant.id = assoc$variant.id)
  assoc$chr <- as.numeric(assoc$chr)
  assoc_out <- merge(assoc, ea, by = "variant.id")
  assoc_out_sorted <- assoc_out[order(assoc_out$chr),]

  # Save output files
  output_file <- paste0(diro, "GENESIS_assoc_", phenotype, "_afr_out.txt")
  write.table(assoc_out_sorted, file = output_file, row.names = FALSE, quote = FALSE)

  # Create Manhattan and Q-Q plot
  manhattan_file <- paste0(diro, "GENESIS_", phenotype, "_assoc_afr_manhattan.png")
  png(file = manhattan_file, width = 8, height = 7, units = "in", res = 300)
  manhattan(assoc_out_sorted, chr = "chr", bp = "pos", snp = "variant.id", p = "Score.pval")
  dev.off()

  qq_file <- paste0(diro, "GENESIS_", phenotype, "_assoc_afr_qq.png")
  png(file = qq_file, width = 8, height = 7, units = "in", res = 300)
  qq(assoc_out_sorted$Score.pval)
  dev.off()
}
