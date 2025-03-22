# Load necessary libraries
library(MOFA2)
library(limma)
library(edgeR)
library(ggplot2)
install.packages("ggdendro")  
install.packages("dendextend")  

ls("package:MOFA2")
??MOFA2

# Define paths to the files
gut_path <- "Gut_Metagenomics_CLR_Normalized.csv"
saliva_path <- "Saliva_Metagenomics_CLR_Normalized.csv"
metabolomics_path <- "normalized_metabolomics.csv"
proteomics_path <- "normalized_proteomics.csv"
metadata_path <- "Metadata.csv"

# Load data
gut <- read.csv(gut_path, row.names = 1)
saliva <- read.csv(saliva_path, row.names = 1)
metabolomics <- read.csv(metabolomics_path, row.names = 1)
proteomics <- read.csv(proteomics_path, row.names = 1)
metadata <- read.csv(metadata_path)

# Transpose if necessary (samples in rows, features in columns)
#gut <- t(gut)
#saliva <- t(saliva)
#metabolomics <- t(metabolomics)
#proteomics <- t(proteomics)

# Normalize data (example using log transformation and Z-score)
#normalize_data <- function(data) {
  #log_data <- log1p(data)  # Log transformation
  #scale(log_data)          # Z-score normalization
#}

#gut <- normalize_data(gut)
#saliva <- normalize_data(saliva)
#metabolomics <- normalize_data(metabolomics)


# Find common samples across all datasets and metadata
common_samples <- Reduce(intersect, list(
  colnames(gut), 
  colnames(saliva), 
  colnames(metabolomics), 
  colnames(proteomics), 
  metadata$Patient.ID
))

# Subset each dataset based on truly common samples
aligned_data <- list(
  Gut = gut[, common_samples, drop = FALSE],
  Saliva = saliva[, common_samples, drop = FALSE],
  Metabolomics = metabolomics[, common_samples, drop = FALSE],
  Proteomics = proteomics[, common_samples, drop = FALSE]
)

# Check if samples are properly aligned
lapply(aligned_data, dim)  # Should now have non-zero rows



# Find common samples across all views
common_samples <- Reduce(intersect, lapply(aligned_data, colnames))

# Subset each view to include only common samples, in the same order
aligned_data <- lapply(aligned_data, function(x) x[, common_samples, drop = FALSE])

# Ensure each dataset is a matrix
aligned_data <- lapply(aligned_data, function(x) as.matrix(x))

# Ensure aligned_data is correctly formatted (features as rows, samples as columns)
str(aligned_data)

# Create MOFA object
mofa <- create_mofa(aligned_data)

# Get default model and training options (requires MOFA object)
model_opts <- get_default_model_options(mofa)
train_opts <- get_default_training_options(mofa)

# Prepare the MOFA model
mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)

# Train MOFA+ model
mofa <- run_mofa(mofa, use_basilisk = TRUE)

# Save MOFA model
saveRDS(mofa, file = "MOFA_model.rds")



# Visualize results
#Factor Variance Explained: Shows how much variance each factor explains in different omics layers.
plot_variance_explained(mofa)

#Factor Scores: Visualizes how samples are distributed along the learned factors.
plot_factors(mofa)

#Factor Associations with Metadata: Helps identify relationships between factors and metadata (e.g., ADAS-Cog score).
samples <- samples_names(mofa)
print(samples)  # Ensure these match your metadata Patient.ID column

print(nrow(metadata))  # Should return 73

mofa_samples <- unlist(samples_names(mofa))  # Convert list to vector
print(length(mofa_samples))  # Should now return 73
print(head(mofa_samples))  # Check sample names

metadata <- metadata[metadata$Patient.ID %in% mofa_samples, ]
metadata <- metadata[match(mofa_samples, metadata$Patient.ID), ]  # Ensure correct order



sum(is.na(metadata$ADAS.group))  # Should be 0
print(nrow(metadata))  # Should return 73
print(all(metadata$Patient.ID == mofa_samples))  # Should return TRUE


plot_factors(mofa, color_by = metadata$ADAS.group)

#Top Features per Factor: Shows the most influential features driving each factor.
plot_top_weights(mofa, view = "Proteomics", factor = 4, nfeatures = 20)
plot_top_weights(mofa, view = "Metabolomics", factor = 4, nfeatures = 20)
plot_top_weights(mofa, view = "Gut", factor = 4, nfeatures = 20)
plot_top_weights(mofa, view = "Saliva", factor = 4, nfeatures = 20)

#Heatmap of Feature Weights: Visualizes how features contribute to each factor.
#plot_weights_heatmap(mofa, view = "Proteomics")

#Sample Clustering Based on Factors: Groups samples using factor scores.
#plot_dendrogram(mofa)

# Extract factor scores (assuming single group 'group1')
factors <- get_factors(mofa)$group1

# Compute a distance matrix between samples (transpose so samples are rows)
d <- dist(t(factors))

# Perform hierarchical clustering
hc <- hclust(d, method = "complete")

# Plot the dendrogram
plot(hc, main = "Dendrogram of Samples based on MOFA Factors")


#Contribution of Each Omics Layer to Factors
plot_data_overview(mofa)

# Latent Space Representation
#UMAP/T-SNE Representation of Samples: Projects samples into a low-dimensional space based on factor scores.
plot_dimred(mofa, method = "UMAP")

#Factor Correlation Heatmap
plot_factor_cor(mofa)

#Bar Plot of Variance Explained
plot_variance_explained(mofa, x = "view", y = "factor")



factors <- get_factors(mofa)
print(str(factors))  # Check if it's a list, matrix, or dataframe
factors <- factors$group1  # Extract the matrix
print(dim(factors))  # Should return (73, 15)
print(head(factors))  # Preview values

aov_test <- aov(factors[,4] ~ metadata$ADAS.group)
summary(aov_test)

anova_results <- apply(factors, 2, function(x) summary(aov(x ~ metadata$ADAS.group))[[1]][["Pr(>F)"]][1])
print(anova_results)  # p-values for each factor

pairwise.t.test(factors[,4], metadata$ADAS.group, p.adjust.method = "bonferroni")
pairwise.wilcox.test(factors[,4], metadata$ADAS.group, p.adjust.method = "bonferroni")

ggplot(metadata, aes(x = ADAS.group, y = factors[,7], fill = ADAS.group)) +
  geom_boxplot() +
  labs(title = "MOFA Factor 7 vs ADAS Group", y = "Factor7 Value") +
  theme_minimal()

ggplot(metadata, aes(x = ADAS.group, y = factors[,4], fill = ADAS.group)) +
  geom_boxplot() +
  labs(title = "MOFA Factor 4 vs ADAS Group", y = "Factor4 Value") +
  theme_minimal()


print(mofa)
