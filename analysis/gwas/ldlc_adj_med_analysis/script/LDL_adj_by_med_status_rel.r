# Script outline
# This script will create multiple LDL cholesterol columns and adjust each column based on different relative medication adjustments

setwd('/user/work/ac14629/MRC_network_project/data/UKB/phenotype')

## Relative LDL-lowering medication adjustment

# Read in data
df <- read.table("mrc_phenotypes_sas.txt", header = T)

# Extract IEU ID and ldlc column
df_ldlc <- df[, c("IID","ieu","med_ldlc","ldlc")]

# Duplicate IEU ID column
id_dup <- matrix(rep(df_ldlc$IID, times = 2), ncol = 2, byrow = FALSE)

# Duplicate ldlc column to create nine ldlc columns
ldlc_dup <- matrix(rep(df_ldlc$ldlc, times = 9), ncol = 9, byrow = FALSE)

# Generate new column names
colnames(id_dup) <- c('FID', 'IID')
colnames(ldlc_dup) <- c('ldlc_unadjusted', 'ldlc_10', 'ldlc_20', 'ldlc_30', 'ldlc_40', 'ldlc_50', 'ldlc_60', 'ldlc_70', 'ldlc_80')

# Assign duplicated columns and new names to data frame
new_df <- data.frame(id_dup, df_ldlc$ieu, df_ldlc$med_ldlc, ldlc_dup)

# Remove NAs from FID and ldl columns
new_df <- new_df[!is.na(new_df$FID), ]
new_df <- new_df[!is.na(new_df$ldlc_unadjusted), ]

# Divide ldlc_10 column value by 0.9 if medication equals 1
new_df$ldlc_10[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_10[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.9

# Divide ldlc_20 column value by 0.8 if medication equals 1
new_df$ldlc_20[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_20[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.8

# Divide ldlc_30 column value by 0.7 if medication equals 1
new_df$ldlc_30[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_30[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.7

# Divide ldlc_40 column value by 0.6 if medication equals 1
new_df$ldlc_40[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_40[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.6

# Divide ldlc_50 column value by 0.5 if medication equals 1
new_df$ldlc_50[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_50[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.5

# Divide ldlc_60 column value by 0.4 if medication equals 1
new_df$ldlc_60[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_60[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.4

# Divide ldlc_70 column value by 0.3 if medication equals 1
new_df$ldlc_70[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_70[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.3

# Divide ldlc_80 column value by 0.2 if medication equals 1
new_df$ldlc_80[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_80[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] / 0.2

# Get the columns with names starting with "ldlc_"
ldlc_cols <- grep("^ldlc_", names(new_df), value = TRUE)

# Apply rank-based inverse normal transformation to ldlc_ columns
for (col in ldlc_cols) {
  ranked_values <- rank(new_df[[col]], na.last = "keep")
  transformed_values <- qnorm((ranked_values - 0.5) / sum(!is.na(ranked_values)))
  new_df[[col]] <- transformed_values
}

# Save merged file
write.table(new_df, "mrc_ldlc_adj_rel_med_sas.txt", sep = "\t", row.names = FALSE)
