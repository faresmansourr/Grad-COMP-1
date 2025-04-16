library(ReactomePA)
library(MetaboAnalystR)
library(MOFA2)            # MOFA analysis
library(clusterProfiler)   # KEGG & GO enrichment
library(org.Hs.eg.db)      # Human gene annotation
library(ReactomePA)        # Reactome pathway enrichment
library(enrichplot)        # Pathway visualization
library(ggplot2)
library(dplyr)
library(readxl)
library(taxize)


# Load the trained MOFA model
pathway_mofa <- readRDS("MOFA_model_scaled.rds")

# Extract factors from MOFA
factors <- get_factors(pathway_mofa)$group1  # Extract factor scores for the first group

# Specify factors of interest
selected_factors <- c("Factor4", "Factor6")

# Get feature loadings (weights) per omics layer
loadings <- get_weights(pathway_mofa)

# Extract features for the selected factors
selected_features <- list()
for (factor in selected_factors) {
  for (view in names(loadings)) {
    df <- as.data.frame(loadings[[view]][, factor, drop=FALSE])
    colnames(df) <- c("Weight")
    df <- df[order(-abs(df$Weight)), , drop=FALSE]  # Rank by absolute weight
    selected_features[[paste(view, factor, sep="_")]] <- df
  }
}

#####################
# Extract top 100 proteins from Factor4 and Factor6
proteomics_features_Factor4 <- rownames(selected_features[["Proteomics_Factor4"]])[1:100]
proteomics_features_Factor6 <- rownames(selected_features[["Proteomics_Factor6"]])[1:100]
proteomics_features <- unique(c(proteomics_features_Factor4, proteomics_features_Factor6))
proteomics_features <- as.character(proteomics_features)
write.csv(proteomics_features, file = 'top_proteins_Factor4_Factor6.csv')

head(proteomics_features)
proteomics_features_clean <- sub("\\..*", "", proteomics_features)
write.csv(proteomics_features_clean, file = 'top_proteins_Factor4_Factor6_clean.csv')

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
dotplot(kegg_enrich, showCategory=15) + ggtitle("KEGG Pathway Enrichment - Proteomics Factors 4 & 6")
dotplot(reactome_enrich, showCategory=15) + ggtitle("Reactome Pathway Enrichment - Proteomics Factors 4 & 6")
dotplot(go_enrich, showCategory=15) + ggtitle("GO Pathway Enrichment - Proteomics Factors 4 & 6")

# Save results
write.csv(kegg_enrich@result, "KEGG_Proteomics_Factor4_Factor6.csv")
write.csv(reactome_enrich@result, "Reactome_Proteomics_Factor4_Factor6.csv")
write.csv(go_enrich@result, "GO_BP_Proteomics_Factor4_Factor6.csv")
##########################

###########################
metabolomics_features_Factor4 <- rownames(selected_features[["Metabolomics_Factor4"]])[1:100]
metabolomics_features_Factor6 <- rownames(selected_features[["Metabolomics_Factor6"]])[1:100]
metabolomics_features <- unique(c(metabolomics_features_Factor4, metabolomics_features_Factor6))
metabolomics_features <- as.character(metabolomics_features)
write.csv(metabolomics_features, file = 'top_metabolites_Factor4_Factor6.csv')
#########################

#########################
gut_features_Factor4 <- rownames(selected_features[["Gut_Factor4"]])[1:100]
gut_features_Factor6 <- rownames(selected_features[["Gut_Factor6"]])[1:100]
gut_features <- unique(c(gut_features_Factor4, gut_features_Factor6))
gut_features <- as.character(gut_features)
write.csv(gut_features, file = 'top_gut_Factor4_Factor6.csv')
#########################

#########################
saliva_features_Factor4 <- rownames(selected_features[["Saliva_Factor4"]])[1:100]
saliva_features_Factor6 <- rownames(selected_features[["Saliva_Factor6"]])[1:100]
saliva_features <- unique(c(saliva_features_Factor4, saliva_features_Factor6))
saliva_features <- as.character(saliva_features)
write.csv(saliva_features, file = 'top_saliva_Factor4_Factor6.csv')
#########################