# Load necessary libraries
library(tidyverse)
library(compositions)  # For CLR transformation
library(ggplot2)
library(limma)

# Define paths to the files
gut_raw_path <- "Gut_Metagenomics.csv"
saliva_raw_path <- "Saliva_Metagenomics.csv"
metabolomics_raw_path <- "Metabolomics.csv"
proteomics_raw_path <- "Proteomics.csv"

# Load data
gut_raw <- read.csv(gut_raw_path, row.names = 1)
saliva_raw <- read.csv(saliva_raw_path, row.names = 1)
metabolomics_raw <- read.csv(metabolomics_raw_path, row.names = 1)
proteomics_raw <- read.csv(proteomics_raw_path, row.names = 1)

gut_raw = t(gut_raw)
saliva_raw = t(saliva_raw)
metabolomics_raw = t(metabolomics_raw)
proteomics_raw = t(proteomics_raw)

# Check if log transformation is needed (should be roughly normal)
hist(as.numeric(proteomics_raw), breaks = 50, main = "Raw Proteomics Data Distribution")

# Boxplot to check batch/sample effects
boxplot(proteomics_raw, main = "Raw Proteomics Data Boxplot", las = 2)

library(preprocessCore)
normalized_proteomics <- normalize.quantiles(as.matrix(proteomics_raw))
boxplot(normalized_proteomics, main="Normalized Proteomics Data Boxplot")
hist(as.numeric(normalized_proteomics), breaks = 50, main = "Normalized Proteomics Data Distribution")

feature_names <- rownames(proteomics_raw)  # Save feature names (proteins/metabolites)
sample_names <- colnames(proteomics_raw)   # Save sample names (patients)
# Restore row and column names
rownames(normalized_proteomics) <- feature_names
colnames(normalized_proteomics) <- sample_names

# Convert back to data frame if needed
normalized_proteomics <- as.data.frame(normalized_proteomics)
# Save the normalized data
write.csv(normalized_proteomics, "normalized_proteomics.csv", row.names = TRUE)
###########################################
# Load necessary libraries

# Check if log transformation is needed (should be roughly normal)
hist(as.numeric(metabolomics_raw), breaks = 50, main = "Raw Metabolomics Data Distribution")

# Boxplot to check batch/sample effects
boxplot(metabolomics_raw, main = "Raw Metabolomics Data Boxplot", las = 2)



# Convert to matrix for normalization
metabolomics_matrix <- as.matrix(metabolomics_raw)

# Log2 transform if needed (check distribution first)
log_metabolomics <- log2(metabolomics_matrix + 1)  # Adding 1 to avoid log(0) issues

# Apply quantile normalization
normalized_metabolomics <- normalize.quantiles(log_metabolomics)

# Restore row and column names
rownames(normalized_metabolomics) <- rownames(metabolomics_matrix)
colnames(normalized_metabolomics) <- colnames(metabolomics_matrix)

# Convert back to data frame if needed
normalized_metabolomics_df <- as.data.frame(normalized_metabolomics)

# Save the normalized data
write.csv(normalized_metabolomics_df, "normalized_metabolomics.csv", row.names = TRUE)

# Boxplot for visualization
boxplot(normalized_metabolomics_df, main = "Normalized Metabolomics Data", las = 2)

###################################
str(gut_raw)
summary(gut_raw)
colSums(gut_raw)
hist(as.numeric(gut_raw), breaks = 50, main = "Raw Gut Metagenomics Data Distribution", xlab = "Values")
# Boxplot to check batch/sample effects
boxplot(gut_raw, main = "Raw Gut Data Boxplot", las = 2)



# Convert to a numeric matrix (assuming taxa are in rows, samples in columns)
data_matrix <- as.matrix(gut_raw)

# Remove taxa with too many zeros (optional, adjust threshold as needed)
zero_threshold <- 0.9  # Remove taxa with >90% zero values
keep_taxa <- rowMeans(data_matrix == 0) < zero_threshold
data_filtered <- data_matrix[keep_taxa, ]

# CLR Transformation (Preferred for MOFA+)
data_clr <- t(apply(data_filtered, 1, function(x) {
  x[x == 0] <- NA  # Avoid log(0) issues
  clr_val <- clr(x, na.rm = TRUE)
  return(clr_val)
}))

# Replace any remaining NA values with zero (optional)
data_clr[is.na(data_clr)] <- 0  

# Alternative: Log-transformation after TSS normalization
data_tss <- sweep(data_filtered, 2, colSums(data_filtered), FUN = "/")  # TSS normalization
data_log <- log2(data_tss + 1)  # Log transform

# Save results
write.csv(data_clr, "Gut_Metagenomics_CLR_Normalized.csv")
write.csv(data_log, "Gut_Metagenomics_LogTSS_Normalized.csv")

# Print message
cat("Normalization complete! Saved CLR-normalized and Log-TSS-normalized data.")
# Boxplot for visualization
boxplot(data_log, main = "Normalized Log-TSS Gut Data", las = 2)
boxplot(data_clr, main = "Normalized CLR Gut Data", las = 2)
##################################
str(saliva_raw)
summary(saliva_raw)
colSums(saliva_raw)
hist(as.numeric(saliva_raw), breaks = 50, main = "Raw Saliva Metagenomics Data Distribution", xlab = "Values")
# Boxplot to check batch/sample effects
boxplot(saliva_raw, main = "Raw Saliva Data Boxplot", las = 2)



# Convert to a numeric matrix (assuming taxa are in rows, samples in columns)
s_data_matrix <- as.matrix(saliva_raw)

# Remove taxa with too many zeros (optional, adjust threshold as needed)
zero_threshold <- 0.9  # Remove taxa with >90% zero values
s_keep_taxa <- rowMeans(s_data_matrix == 0) < zero_threshold
s_data_filtered <- s_data_matrix[s_keep_taxa, ]

# CLR Transformation (Preferred for MOFA+)
s_data_clr <- t(apply(s_data_filtered, 1, function(x) {
  x[x == 0] <- NA  # Avoid log(0) issues
  clr_val <- clr(x, na.rm = TRUE)
  return(clr_val)
}))

# Replace any remaining NA values with zero (optional)
s_data_clr[is.na(s_data_clr)] <- 0  

# Alternative: Log-transformation after TSS normalization
s_data_tss <- sweep(s_data_filtered, 2, colSums(s_data_filtered), FUN = "/")  # TSS normalization
s_data_log <- log2(s_data_tss + 1)  # Log transform

# Save results
write.csv(s_data_clr, "Saliva_Metagenomics_CLR_Normalized.csv")
write.csv(s_data_log, "Saliva_Metagenomics_LogTSS_Normalized.csv")

# Print message
cat("Normalization complete! Saved CLR-normalized and Log-TSS-normalized data.")
# Boxplot for visualization
boxplot(s_data_log, main = "Normalized Log-TSS Saliva Data", las = 2)
boxplot(s_data_clr, main = "Normalized CLR Saliva Data", las = 2)
