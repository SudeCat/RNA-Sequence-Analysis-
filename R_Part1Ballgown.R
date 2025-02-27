# Load required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data
pheno_data <- read.csv("comparison.csv")

# Create Ballgown object
bg <- ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Perform differential expression analysis
results_transcripts <- stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes <- stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

# Save results
write.table(results_transcripts, "transcript_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(results_genes, "gene_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Filter low-abundance genes
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Perform DE analysis on filtered data
results_genes_filtered <- stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

# Save filtered results
write.table(results_genes_filtered, "gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Identify significant genes (p-value < 0.05)
sig_genes <- subset(results_genes_filtered, results_genes_filtered$pval < 0.05)
write.table(sig_genes, "gene_results_sig.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Exit the R session
quit(save="no")