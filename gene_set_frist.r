if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("msigdbr")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
# Example gene list and logFC
geneList <- c("TP53", "BRCA1", "BRCA2", "RB1", "PTEN", "MDM2", "CDK2", "CDK4")
logFC <- c(-1.5, 2.4, -2.1, 1.8, -3.2, 2.5, -1.7, 1.9)

# Convert gene symbols to Entrez IDs
entrezIDs <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
names(entrezIDs) <- entrezIDs$ENTREZID
entrezList <- setNames(entrezIDs$ENTREZID, entrezIDs$SYMBOL)

# Rank genes
rankedGenes <- sort(entrezList, decreasing = TRUE)
# Load MSigDB ontology gene sets for humans
msigdb <- msigdbr(species = "Homo sapiens", category = "C5")
gseaResults <- GSEA(geneList = rankedGenes, 
                    exponent = 1,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pAdjustMethod = "BH",
                    TERM2GENE = msigdb[, c("gs_name", "gene_symbol")])
# Summary of results
summary(gseaResults)

# Enrichment plot for the most significant gene set
enrichPlot(gseaResults, showCategory = 5) + ggplot2::theme_minimal()

# Dotplot for visualizing multiple gene sets
dotplot(gseaResults) + ggplot2::theme_minimal()
