library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)
library(enrichplot)

filename <- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)
my_symbols_20 <- c("HSP90B1", "F8", "SIAH2", "GLUL", "SLC12A9",
                   "COL3A1", "GJB3", "MYH9", "ITGB1", "LAMC2",
                   "CERCAM", "SPINT2", "SLC2A1", "LDHA", "CTSC",
                   "G6PD", "LAMA5", "KRT17", "RASD2", "SEZ6L")

hs <- org.Hs.eg.db
entrez_20 <- select(hs, keys = my_symbols_20, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
entrez_20 <- na.omit(entrez_20)  # NA değerleri kaldır

ImmunSigEnrich_20 <- enricher(entrez_20$ENTREZID, TERM2GENE = c6, pvalueCutoff = 0.05)
ImmunSigEnrich_20 <- setReadable(ImmunSigEnrich_20, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(ImmunSigEnrich_20, "ANNOVAR_ImmunePathway.csv")

go_20 <- enrichGO(gene = entrez_20$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont = "ALL", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, readable = TRUE)

write.csv(go_20, "ANNOVAR_GO_Enrichment.csv")

kegg_genes_20 <- bitr(entrez_20$ENTREZID, fromType = "ENTREZID", toType = "KEGG", OrgDb = org.Hs.eg.db)
valid_kegg_genes_20 <- na.omit(kegg_genes_20$KEGG)
kegg_20 <- enrichKEGG(gene = valid_kegg_genes_20, organism = "hsa",

pAdjustMethod = "BH", pvalueCutoff = 0.1)

write.csv(kegg_20, "ANNOVAR_KEGG_Enrichment.csv")
pdf(file = "ANNOVAR_GSEA_Plots.pdf")

barplot(go_20, showCategory = 20)
dotplot(go_20, showCategory = 30)
barplot(kegg_20, showCategory = 20)
dotplot(kegg_20, showCategory = 30)

pdf("Annovar_Cytoscape_Network.pdf") 
cnetplot(go_annovar, showCategory = 5, foldChange = annovar_entrez)

# Combine both gene lists
combined_entrez <- unique(c(htseq_entrez, annovar_entrez))

# Perform GO Enrichment Analysis for combined genes
go_combined <- enrichGO(gene = combined_entrez, OrgDb = org.Hs.eg.db, ont = "ALL",
pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# Save Combined GO Enrichment results
write.csv(go_combined, "Combined_GO_Enrichment.csv")

# Generate Combined Cytoscape Network Plot
pdf("Combined_Cytoscape_Network.pdf")
cnetplot(go_combined, showCategory = 5, foldChange = combined_entrez)
dev.off()	
