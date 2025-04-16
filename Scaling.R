# Load necessary libraries
library(data.table)
library(ggplot2)

# Function to apply z-score scaling
df_zscore <- function(df) {
  return(scale(df))  # Centers and scales each feature
}

# Load normalized datasets
gut_metagenomics_norm <- fread("Gut_Metagenomics_CLR_Normalized.csv", header = TRUE)
saliva_metagenomics_norm <- fread("Saliva_Metagenomics_CLR_Normalized.csv", header = TRUE)
proteomics_norm <- fread("normalized_proteomics.csv", header = TRUE)
metabolomics_norm <- fread("normalized_metabolomics.csv", header = TRUE)

# Ensure first column has a proper name
colnames(gut_metagenomics_norm)[1] <- "Feature"
colnames(saliva_metagenomics_norm)[1] <- "Feature"
colnames(proteomics_norm)[1] <- "Feature"
colnames(metabolomics_norm)[1] <- "Feature"

# Convert to numeric matrices (excluding Feature column)
gut_metagenomics_scaled <- as.data.frame(df_zscore(as.matrix(gut_metagenomics_norm[,-1])))
saliva_metagenomics_scaled <- as.data.frame(df_zscore(as.matrix(saliva_metagenomics_norm[,-1])))
proteomics_scaled <- as.data.frame(df_zscore(as.matrix(proteomics_norm[,-1])))
metabolomics_scaled <- as.data.frame(df_zscore(as.matrix(metabolomics_norm[,-1])))

# Add back Feature column
gut_metagenomics_scaled <- cbind(Feature = gut_metagenomics_norm$Feature, gut_metagenomics_scaled)
saliva_metagenomics_scaled <- cbind(Feature = saliva_metagenomics_norm$Feature, saliva_metagenomics_scaled)
proteomics_scaled <- cbind(Feature = proteomics_norm$Feature, proteomics_scaled)
metabolomics_scaled <- cbind(Feature = metabolomics_norm$Feature, metabolomics_scaled)

# Save scaled data
fwrite(gut_metagenomics_scaled, "Gut_Metagenomics_Scaled.csv")
fwrite(saliva_metagenomics_scaled, "Saliva_Metagenomics_Scaled.csv")
fwrite(proteomics_scaled, "Proteomics_Scaled.csv")
fwrite(metabolomics_scaled, "Metabolomics_Scaled.csv")

#Visualization
boxplot(gut_metagenomics_scaled[, -1], 
        main = "Gut Metagenomics Scaled Boxplot", 
        las = 2)

boxplot(saliva_metagenomics_scaled[, -1], 
        main = "Saliva Metagenomics Scaled Boxplot", 
        las = 2)

boxplot(proteomics_scaled[, -1], 
        main = "Proteomics Scaled Boxplot", 
        las = 2)

boxplot(metabolomics_scaled[, -1], 
        main = "Metabolomics Scaled Boxplot", 
        las = 2)

