### Load Required Libraries ###
library(gplots) 
library(RColorBrewer) 

### Define Global Parameters ###
min_nonzero <- 1e-6 # Avoid log2(0) issues
data_columns <- grep("FPKM", colnames(gene_expression)) # Identify FPKM columns dynamically
short_names <- gsub("FPKM\\.|_", "", colnames(gene_expression)[data_columns]) # Shorten column names
data_colors <- rainbow(length(data_columns)) # Assign colors to samples

### Plot #1 - Number of Transcripts Per Gene ###
counts <- table(transcript_gene_table[,"g_id"])
c_one <- length(which(counts == 1))
c_more_than_one <- length(which(counts > 1))
c_max <- max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene",
     main="Distribution of transcript count per gene")
legend("topright", legend=c(paste("Genes with one transcript =", c_one), 
                            paste("Genes with more than one transcript =", c_more_than_one), 
                            paste("Max transcripts for single gene =", c_max)), lty=NULL)

### Plot #2 - Transcript Length Distribution ###
full_table <- texpr(bg, 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", 
main="Distribution of transcript lengths", col="steelblue")

### Plot #3 - Distribution of FPKM Values Across Libraries ###
boxplot(log2(gene_expression[, data_columns] + min_nonzero), 
        col=data_colors, names=short_names, las=2, ylab="log2(FPKM)",
        main="Distribution of FPKMs for all libraries")

### Plot #4 - Pairwise Correlation Between Two Replicates ###
x <- gene_expression[, "FPKM.UHR_Rep57"]
y <- gene_expression[, "FPKM.UHR_Rep58"]
plot(log2(x + min_nonzero), log2(y + min_nonzero), pch=16, col="blue", cex=0.25,
     xlab="FPKM (UHR, Replicate 57)", ylab="FPKM (UHR, Replicate 58)",
     main="Comparison of expression values for a pair of replicates")
abline(a=0, b=1)
rs <- cor(x, y, use="complete.obs")^2
legend("topleft", paste("R squared =", round(rs, 3)), lwd=1, col="black")

### Plot #5 - Density Scatter Plot ###
colors <- colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(log2(x + min_nonzero), log2(y + min_nonzero), 
              xlab="FPKM (UHR, Replicate 57)", ylab="FPKM (UHR, Replicate 58)",
              main="Density Scatter Plot: UHR Replicates 57 and 58",
              colramp=colors, nbin=200)
abline(a=0, b=1)

### Plot #6 - Multiple Replicate Comparisons ###
plotCor <- function(lib1, lib2, name, color) {
  x <- gene_expression[, lib1]
  y <- gene_expression[, lib2]
  zero_count <- length(which(x == 0)) + length(which(y == 0))
 plot(log2(x + min_nonzero), log2(y + min_nonzero), pch=16, col=color, cex=0.25,
       xlab=lib1, ylab=lib2, main=name)abline(a=0, b=1, col="red")
rs <- cor(x, y, method="pearson")^2
legend("topleft", c(paste("R² =", round(rs, 3)), paste("Zero count =", zero_count)), col="black")
}
par(mfrow=c(1,2))
plotCor("FPKM.UHR_Rep57", "FPKM.HBR_Rep63", "UHR_57 vs HBR_63", "tomato2")
plotCor("FPKM.UHR_Rep58", "FPKM.HBR_Rep64", "UHR_58 vs HBR_64", "royalblue2")

### Plot #7 - Multidimensional Scaling (MDS) ###
d <- 1 - cor(gene_expression[, data_columns], use="pairwise.complete.obs", method="pearson")
mds <- cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", main="MDS Distance Plot (All Non-Zero Genes)")
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

### Plot #8 - Distribution of Differential Expression ###
sig <- which(results_genes$pval < 0.05)
results_genes[, "de"] <- log2(results_genes[, "fc"])
hist(results_genes[sig, "de"], breaks=50, col="seagreen",
     xlab="log2(Fold Change) UHR vs HBR", main="Differential Expression Distribution")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)

### Plot #9 - Heatmap for Differentially Expressed Genes ###
mydist <- function(c) { dist(c, method="euclidean") }
myclust <- function(c) { hclust(c, method="average") }
sig_genes <- results_genes[sig, "id"]
sig_gene_names <- results_genes[sig, "gene_name"]
data <- log2(as.matrix(gene_expression[sig_genes, data_columns]) + 1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm=TRUE, scale="none",
          dendrogram="both", margins=c(6,7), Rowv=TRUE, Colv=TRUE,
          key=TRUE, trace="none", main="Heatmap of Significant DE Transcripts",
          cexRow=0.3, cexCol=1, labRow=sig_gene_names, col=rev(heat.colors(75)))

### Save Differentially Expressed Genes ###
sigpi <- which(results_genes$pval < 0.05)
sigp <- results_genes[sigpi,]
sigde <- which(abs(sigp[, "de"]) >= 2)
sig_tn_de <- sigp[sigde,]
o <- order(sig_tn_de[, "qval"], -abs(sig_tn_de[, "de"]), decreasing=FALSE)
output <- sig_tn_de[o, c("gene_name", "id", "fc", "pval", "qval", "de")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Exit the R session
quit(save="no")