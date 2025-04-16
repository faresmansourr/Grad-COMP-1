library(readxl)
library(taxize)
library(rentrez)
library(ggplot2)

# === Step 1: Load Taxa Data from Excel ===
file_path <- "Top_Taxa_Factors4,6.xlsx"
taxa_data <- readxl::read_excel(file_path)

taxa_list <- na.omit(taxa_data$Top_taxa_all)
taxa_list <- unique(taxa_list)

# === Step 2: Map Taxa to NCBI Taxonomy IDs ===
get_taxid_safe <- function(taxon) {
  tryCatch({
    res <- taxize::get_uid(taxon, messages = FALSE)
    if (length(res[[1]]) > 0) return(res[[1]][1]) else return(NA)
  }, error = function(e) NA)
}

cat("Mapping taxa to NCBI TaxIDs...\n")
taxid_results <- data.frame(
  Taxon = taxa_list,
  TaxID = sapply(taxa_list, get_taxid_safe),
  stringsAsFactors = FALSE
)

write.csv(taxid_results, "MOFA_Taxa_to_TaxIDs.csv", row.names = FALSE)
cat("Taxon-to-TaxID mapping complete. Results saved to 'MOFA_Taxa_to_TaxIDs.csv'.\n")

# === Step 3: Fetch Genome Assembly Accessions from NCBI ===
get_assembly_accession <- function(taxid) {
  if (is.na(taxid)) return(NA)
  tryCatch({
    search <- entrez_search(db = "assembly", term = paste0(taxid, "[TaxID]"), retmax = 1)
    if (length(search$ids) > 0) {
      summary <- entrez_summary(db = "assembly", id = search$ids[[1]])
      return(summary$assemblyaccession)
    } else return(NA)
  }, error = function(e) NA)
}

cat("Fetching genome assembly accessions...\n")
taxid_results$Assembly_Accession <- sapply(taxid_results$TaxID, get_assembly_accession)

write.csv(taxid_results, "MOFA_Taxa_to_Assemblies.csv", row.names = FALSE)
cat("Assembly accessions added. Results saved to 'MOFA_Taxa_to_Assemblies.csv'.\n")

# === Step 4: Summary Plot (Optional) ===
ggplot(taxid_results, aes(x = !is.na(Assembly_Accession))) +
  geom_bar(fill = "steelblue") +
  labs(title = "Number of Taxa with Genome Assemblies", x = "Has Assembly Accession", y = "Count") +
  theme_minimal()

ggsave("Taxa_Genome_Assembly_Summary.png")


