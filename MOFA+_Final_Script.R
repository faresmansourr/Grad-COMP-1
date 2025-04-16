# Load necessary libraries
library(MOFA2)
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(data.table)
library(GGally)  # Load the package
library(psych)

ls("package:MOFA2")
??MOFA2

# Define paths to the files
gut_path <- "Gut_Metagenomics_Scaled.csv"
saliva_path <- "Saliva_Metagenomics_Scaled.csv"
metabolomics_path <- "Metabolomics_Scaled.csv"
proteomics_path <- "Proteomics_Scaled.csv"
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
saveRDS(mofa, file = "MOFA_model_scaled.rds")



# Visualize results
#Factor Variance Explained: Shows how much variance each factor explains in different omics layers.
plot_variance_explained(mofa, max_r2 = 5)

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

plot_top_weights(mofa, view = "Proteomics", factor = 6, nfeatures = 20)
plot_top_weights(mofa, view = "Metabolomics", factor = 6, nfeatures = 20)
plot_top_weights(mofa, view = "Gut", factor = 6, nfeatures = 20)
plot_top_weights(mofa, view = "Saliva", factor = 6, nfeatures = 20)
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


# Load MOFA model
mofa_model <- readRDS("MOFA_model_scaled.rds")  
metadata <- read.csv("Metadata.csv", row.names = 1)

# Extract factor values
factors <- get_factors(mofa_model, factors = "all")$group1  # Assuming one group; adjust if needed

# Ensure matching sample names
common_samples2 <- intersect(rownames(factors), rownames(metadata))
factors <- factors[common_samples2, , drop = FALSE]
metadata <- metadata[common_samples2, , drop = FALSE]
metadata_numeric <- metadata
metadata_numeric$ADAS.group <- as.numeric(as.factor(metadata$ADAS.group))  # Convert ADAS groups (Low, Medium, High) to 1,2,3

# Compute Spearman correlation between factors and metadata
cor_matrix <- cor(factors, metadata_numeric, method = "spearman", use = "pairwise.complete.obs")

cor.mtest <- function(mat1, mat2, method = "spearman") {
  p.mat <- matrix(NA, ncol = ncol(mat2), nrow = ncol(mat1))
  for (i in 1:ncol(mat1)) {
    for (j in 1:ncol(mat2)) {
      test <- cor.test(mat1[, i], mat2[, j], method = method, exact = FALSE)  # Suppress warning
      p.mat[i, j] <- test$p.value
    }
  }
  colnames(p.mat) <- colnames(mat2)
  rownames(p.mat) <- colnames(mat1)
  return(p.mat)
}



factors <- factors[, apply(factors, 2, function(x) length(unique(x)) > 1)]
metadata_numeric <- metadata_numeric[, apply(metadata_numeric, 2, function(x) length(unique(x)) > 1)]

metadata_numeric[] <- lapply(metadata_numeric, as.numeric)

table(metadata_numeric$ADAS.group)

metadata_numeric$ADAS.group <- as.numeric(as.factor(metadata_numeric$ADAS.group))
colnames(metadata_numeric)


# Generate p-values
p_values <- cor.mtest(factors, metadata_numeric, method = "spearman")

# Visualize correlation matrix
ggcorrplot(cor_matrix, 
           type = "lower", 
           lab = TRUE, 
           p.mat = p_values,  # ✅ Use directly
           insig = "blank", 
           title = "MOFA Factors vs. Metadata Correlation", 
           colors = c("blue", "white", "red"))  


# Ensure ADAS.group is a factor
metadata_numeric$ADAS.group <- as.factor(metadata_numeric$ADAS.group)

# Compute correlation matrix (excluding ADAS.group)
numeric_vars <- metadata_numeric[, !colnames(metadata_numeric) %in% "ADAS.group"]
cor_matrix <- cor(factors, numeric_vars, method = "spearman")

# Compute correlation p-values
p_values <- cor.mtest(factors, numeric_vars, method = "spearman")

# Compute ANOVA p-values for ADAS.group separately
anova_results <- apply(factors, 2, function(f) summary(aov(f ~ metadata_numeric$ADAS.group))[[1]][["Pr(>F)"]][1])

# Convert ANOVA results to a dataframe
anova_df <- data.frame(Factor = names(anova_results), p_value = anova_results)

# Add ANOVA results to correlation matrix (optional)
anova_df$p_value[anova_df$p_value < 0.05] <- -0.5  # Assign a negative value for visualization
anova_matrix <- matrix(NA, nrow = ncol(factors), ncol = 1)
rownames(anova_matrix) <- colnames(factors)
colnames(anova_matrix) <- "ADAS.group"
anova_matrix[, 1] <- anova_df$p_value

# Combine correlation and ANOVA results
final_matrix <- cbind(anova_matrix, cor_matrix)

# Visualization
ggcorrplot(final_matrix, 
           type = "lower", 
           lab = TRUE, 
           insig = "blank", 
           title = "MOFA Factors vs. Metadata Correlation (with ANOVA for ADAS.group)", 
           colors = c("blue", "white", "red"))








# Plot variance explained by factors
plot_variance_explained(mofa, x = "view", y = "factor") +
  ggtitle("Variance Explained by MOFA Factors") +
  theme_minimal()

# Extract MOFA factor scores
factors <- get_factors(mofa)$group1  # Assuming single group 'group1'

# Perform ANOVA to check association with ADAS groups
anova_results <- apply(factors, 2, function(x) {
  summary(aov(x ~ metadata$ADAS.group))[[1]][["Pr(>F)"]][1]
})

# Print p-values for each factor
print(anova_results)

# Adjust for multiple comparisons using Bonferroni correction
p_adjusted <- p.adjust(anova_results, method = "bonferroni")
print(p_adjusted)

# Compute variance of factors
factor_variance <- apply(factors, 2, var)

# Define high and low variance threshold (e.g., median split)
threshold <- median(factor_variance)

# Create dataframe for plotting
factor_df <- data.frame(Factor = colnames(factors), Variance = factor_variance)
factor_df$VarianceGroup <- ifelse(factor_variance > threshold, "High Variance", "Low Variance")

# Plot factor variances
ggplot(factor_df, aes(x = reorder(Factor, Variance), y = Variance, fill = VarianceGroup)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Variance of MOFA Factors", x = "Factor", y = "Variance") +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal()



# Plot MOFA factors associated with ADAS-Cog groups
plot_factor(mofa_model, 
            factors = 1:6, 
            color_by = "ADAS.group", 
            shape_by = "ADAS.group") + 
  ggtitle("MOFA Factors Associated with ADAS-Cog Groups")






# Check if metadata sample names match the MOFA model
print(head(samples_names(mofa_model)))  # Get sample names in MOFA model
print(head(metadata))  # Ensure metadata contains correct sample names

# Assign metadata to MOFA model
samples_metadata(mofa_model) <- metadata

colnames(metadata)  # List column names in metadata
head(metadata)  # Preview the first few rows
metadata$sample <- rownames(metadata)  # Move row names to a new column
rownames(metadata) <- NULL  # Reset row names

# Beeswarm plot for individual factors
p <- plot_factor(mofa_model,
                 factors = 1:6,
                 color_by = "ADAS.group",
                 dot_size = 3,
                 dodge = TRUE,
                 legend = FALSE,
                 add_violin = TRUE,
                 violin_alpha = 0.5) 

# Customizing colors
p <- p + scale_color_manual(values = c("A" = "green", "B" = "red")) +
  scale_fill_manual(values = c("A" = "green", "B" = "red"))

# Print plot
print(p)

# Scatter plot of factor combinations
plot_factors(mofa_model, 
             factors = 1:6, 
             color_by = "ADAS.group") + 
  ggtitle("Scatter Plot of MOFA Factor Combinations")

plot_factor(mofa_model, factors = 1:3, color_by = "ADAS.group") + 
  ggtitle("MOFA Factors by ADAS-Cog Groups")




plot_factor(mofa, 
            factors = 6, 
            color_by = "Factor6"
)
