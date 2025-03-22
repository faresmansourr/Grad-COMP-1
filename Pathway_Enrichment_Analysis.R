BiocManager::install("ReactomePA")
if (!requireNamespace("reactome.db", quietly = TRUE)) {
  BiocManager::install("reactome.db")
}
BiocManager::install(c("RBGL", "impute", "pcaMethods", "MSnbase", "siggenes"))
Sys.setenv(PATH = paste("C:/rtools44/usr/bin;", Sys.getenv("PATH"), sep=""))
install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)

# Try installing MetaboAnalystR again
devtools::install_github("xia-lab/MetaboAnalystR")
system("gfortran --version")

library(ReactomePA)
library(MetaboAnalystR)
library(MOFA2)            # MOFA analysis
library(clusterProfiler)   # KEGG & GO enrichment
library(org.Hs.eg.db)      # Human gene annotation
library(ReactomePA)        # Reactome pathway enrichment
library(enrichplot)        # Pathway visualization
library(ggplot2)
library(dplyr)

ls("package:MetaboAnalystR")

dim(factors)  # Check number of rows in factors
dim(metadata)  # Check number of rows in metadata
length(metadata$ADAS.group)  # Check length of ADAS.group

# Load the trained MOFA model
pathway_mofa <- readRDS("MOFA_model.rds")

# Extract factors from MOFA
factors <- get_factors(mofa)$group1  # Extract factor scores for the first group

# Perform ANOVA to find factors significantly associated with ADAS groups
anova_results <- apply(factors, 2, function(x) summary(aov(x ~ metadata$ADAS.group))[[1]][["Pr(>F)"]][1])

# Select significant factors (p < 0.05)
sig_factors <- names(anova_results[anova_results < 0.05])

# Print significant factors
print(sig_factors)

# Get feature loadings (weights) per omics layer
loadings <- get_weights(mofa)

# Extract features for each significant factor
selected_features <- list()
for (factor in sig_factors) {
  for (view in names(loadings)) {
    df <- as.data.frame(loadings[[view]][, factor, drop=FALSE])
    colnames(df) <- c("Weight")
    df <- df[order(-abs(df$Weight)), , drop=FALSE]  # Rank by absolute weight
    selected_features[[paste(view, factor, sep="_")]] <- df
  }
}

#####################
# Extract top 100 proteins from a significant factor
proteomics_features <- rownames(selected_features[["Proteomics_Factor4"]])[1:100]
proteomics_features <- as.character(proteomics_features)

head(proteomics_features)
proteomics_features_clean <- sub("\\..*", "", proteomics_features)

valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sum(proteomics_features_clean %in% valid_symbols)  # Count valid symbols

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(proteomics_features_clean, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Run KEGG enrichment
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", pAdjustMethod = "BH")

# Run Reactome enrichment
reactome_enrich <- enrichPathway(gene = entrez_ids$ENTREZID, organism = "human", pAdjustMethod = "BH")

# Run GO enrichment (Biological Process)
go_enrich <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")

# Visualize KEGG pathways
dotplot(kegg_enrich, showCategory=15) + ggtitle("KEGG Pathway Enrichment - Proteomics")
dotplot(reactome_enrich, showCategory=15) + ggtitle("Reactome Pathway Enrichment - Proteomics")
dotplot(go_enrich, showCategory=15) + ggtitle("GO Pathway Enrichment - Proteomics")
# Save results
write.csv(kegg_enrich@result, "KEGG_Proteomics_Factor.csv")
write.csv(reactome_enrich@result, "Reactome_Proteomics_Factor.csv")
write.csv(go_enrich@result, "GO_BP_Proteomics_Factor.csv")





#########################
# Extract top 100 metabolites from a significant factor
#metabolomics_features <- rownames(selected_features[["Metabolomics_Factor4"]])[1:125]
#metabolomics_features2 = as.matrix(metabolomics_features)
#write.csv(metabolomics_features2, file = 'top125metabolites.csv')
# Create a new MetaboAnalyst object for metabolomics analysis
# Initialize MetaboAnalyst object
#mSetObj <- InitDataObjects("metpa", "pathway", FALSE)

# Match metabolites (use your feature list)
#mSetObj <- MatchMetabolites(mSetObj, metabolomics_features)
#mSetObj <- MetaboliteMappingExact(mSetObj, q.type = "name")

#lsf.str("package:MetaboAnalystR")  # If using MetaboAnalystR

# Retrieve mapped KEGG IDs
#kegg_mapping <- mSetObj$dataSet$type

# View results
#head(kegg_mapping)

# Convert to KEGG IDs if available
#kegg_metabolites <- read.csv("Metabolite_KEGG_Mapping.csv")  # Ensure it has "Metabolite" and "KEGG_ID"
#metabolomics_kegg <- kegg_metabolites$KEGG_ID[kegg_metabolites$Metabolite %in% metabolomics_features]

# Run pathway enrichment
#kegg_enrich_metabolomics <- enrichKEGG(gene = metabolomics_kegg, organism = "hsa", pAdjustMethod = "BH")

# Visualize
#dotplot(kegg_enrich_metabolomics, showCategory=15) + ggtitle("KEGG Pathway Enrichment - Metabolomics")

# Save results
#write.csv(kegg_enrich_metabolomics@result, "KEGG_Metabolomics_Factor4.csv")





##################
# Extract top taxa from a significant factor
#metagenomics_features <- rownames(selected_features[["Gut_Factor4"]])[1:100]
#metagenomics_features
#write.csv(metagenomics_features, file = 'top100taxa.csv')
# Load predicted KO abundances from PICRUSt2 or Tax4Fun
#ko_table <- read.csv("PICRUSt2_KO.csv", row.names=1)

# Run KEGG enrichment
#kegg_enrich_metagenomics <- enrichKEGG(gene = rownames(ko_table), organism = "hsa", pAdjustMethod = "BH")

# Visualize
#dotplot(kegg_enrich_metagenomics, showCategory=15) + ggtitle("KEGG Pathway Enrichment - Metagenomics")

# Save results
#write.csv(kegg_enrich_metagenomics@result, "KEGG_Metagenomics_Factor6.csv")



###############
#s_metagenomics_features <- rownames(selected_features[["Saliva_Factor4"]])[1:100]
#s_metagenomics_features
#write.csv(s_metagenomics_features, file = 'top100taxa-saliva.csv')


