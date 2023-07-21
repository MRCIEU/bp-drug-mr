# Script outline
# This script will create multiple LDL cholesterol columns and adjust each column based on different absolute medication adjustments

setwd('/user/work/ac14629/MRC_network_project/data/UKB/phenotype')

## Relative LDL-lowering medication adjustment

# Read in data
df <- read.table("mrc_phenotypes_eas.txt", header = T)

# Extract IEU ID and ldlc column
df_ldlc <- df[, c("IID","med_ldlc","ldlc")]

# Duplicate IEU ID column
id_dup <- matrix(rep(df_ldlc$IID, times = 2), ncol = 2, byrow = FALSE)

# Duplicate ldlc column to create nine ldlc columns
ldlc_dup <- matrix(rep(df_ldlc$ldlc, times = 10), ncol = 10, byrow = FALSE)

# Generate new column names
colnames(id_dup) <- c('FID', 'IID')
colnames(ldlc_dup) <- c('ldlc_0.3', 'ldlc_0.5', 'ldlc_0.7', 'ldlc_0.9', 'ldlc_1.0', 'ldlc_1.2', 'ldlc_1.4', 'ldlc_1.6', 'ldlc_1.8', 'ldlc_2.0')

# Assign duplicated columns and new names to data frame
new_df <- data.frame(id_dup, df_ldlc$med_ldlc, ldlc_dup)
new_df <- new_df[!is.na(new_df$ldlc_0.3), ]

# Remove NAs from FID and ldl columns
new_df <- new_df[!is.na(new_df$FID), ]

# Add 0.3 to ldlc_0.3 column value if medication equals 1
new_df$ldlc_0.3[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_0.3[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 0.3

# Add 0.5 to ldlc_0.5 column value if medication equals 1
new_df$ldlc_0.5[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_0.5[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 0.5

# Add 0.7 to ldlc_0.7 column value if medication equals 1
new_df$ldlc_0.7[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_0.7[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 0.7

# Add 0.9 to ldlc_0.9 column value if medication equals 1
new_df$ldlc_0.9[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_0.9[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 0.9

# Add 1.0 to ldlc_1.0 column value if medication equals 1
new_df$ldlc_1.0[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_1.0[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 1

# Add 1.2 to ldlc_1.2 column value if medication equals 1
new_df$ldlc_1.2[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_1.2[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 1.2

# Add 1.4 to ldlc_1.4 column value if medication equals 1
new_df$ldlc_1.4[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_1.4[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 1.4

# Add 1.6 to ldlc_1.6 column value if medication equals 1
new_df$ldlc_1.6[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_1.6[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 1.6

# Add 1.8 to ldlc_1.8 column value if medication equals 1
new_df$ldlc_1.8[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_1.8[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 1.8

# Add 2.0 to ldlc_2.0 column value if medication equals 1
new_df$ldlc_2.0[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] <- new_df$ldlc_2.0[new_df$df_ldlc.med_ldlc == 1 & !is.na(new_df$df_ldlc.med_ldlc)] + 2

# Get the columns with names starting with "ldlc_"
ldlc_cols <- grep("^ldlc_", names(new_df), value = TRUE)

# Apply rank-based inverse normal transformation to ldlc_ columns
for (col in ldlc_cols) {
  ranked_values <- rank(new_df[[col]], na.last = "keep")
  transformed_values <- qnorm((ranked_values - 0.5) / sum(!is.na(ranked_values)))
  new_df[[col]] <- transformed_values
}

# Save merged file
write.table(new_df, "mrc_ldlc_adj_abs_med_eas.txt", sep = "\t", row.names = FALSE)




