# Load required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
setwd("/workspace/rnaseq/de/ballgown/ref_only")

# Load phenotype data
pheno_data = read.csv("UHR_vs_HBR.csv")
print(pheno_data)

# Load the Ballgown object
load('bg.rda')
print(ls())
print(bg)

# Extract FPKM values from the Ballgown object
fpkm = texpr(bg, meas = "FPKM") 
fpkm = log2(fpkm + 1) 
print(tail(fpkm))

# Create a boxplot for all sample FPKM values
boxplot(fpkm, col = as.numeric(pheno_data$type) + 1, las = 2, ylab = 'log2(FPKM+1)')
gene_names = ballgown::geneNames(bg)

# Identify the index of the RASD2 gene in the datasetrasd2_index = which(gene_names == "RASD2")
if (length(rasd2_index) == 0) {
  stop("RASD2 gene is not found in Ballgown object! ")
}
# Get transcript names of RASD2transcript_names = ballgown::transcriptNames(bg)[rasd2_index]
print(paste("RASD2 Gene Transcript ID's: ", paste(transcript_names, collapse = ", ")))
plot(fpkm[rasd2_index, ] ~ pheno_data$type, border = c(2, 3), main = paste("RASD2 Gene: ", paste(transcript_names, collapse = ", ")), pch = 19, xlab = "Type", ylab = "log2(FPKM+1)") 
points(fpkm[rasd2_index, ] ~ jitter(as.numeric(pheno_data$type)), col = as.numeric(pheno_data$type) + 1, pch = 16)

for (sample in unique(pheno_data$ids)) {
plotTranscripts(ballgown::geneIDs(bg)[rasd2_index], bg,main = paste("Transcript Structure for", sample), sample = sample) }

outfile = "/workspace/rnaseq/de/ballgown/ref_only/RASD2_Analysis.pdf"

pdf(file = outfile)

boxplot(fpkm, col = as.numeric(pheno_data$type) + 1, las = 2, ylab = 'log2(FPKM+1)')

dev.off()

quit(save = "no")
