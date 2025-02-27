# Set working directory
setwd("/media/data07/sude/workspace/rnaseq/de/htseq_counts") 

# Read gene mapping mapping 
<- read.table("ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1) 

# Read count matrix 
rawdata <- read.table("gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1) 

# Create class labels 
class <- factor(c(rep("UHR", 6), rep("HBR", 6)))

# Require at least 25% of samples to have count > 25
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 25) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)

# load edgeR
library('edgeR')

# make class labels
class <- factor( c( rep("UHR",3), rep("HBR",3) ))

# Get common gene names
genes=rownames(rawdata)
gene_names=mapping[genes,1]

# Make DGEList object
y <- DGEList(counts=rawdata, genes=genes, group=class)
nrow(y)

# TMM Normalization
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Differential expression test
et <- exactTest(y)

# Print top genes
topTags(et)

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(y)[as.logical(de)]

# Output DE genes
# Matrix of significantly DE genes
mat <- cbind(
 genes,gene_names,
 sprintf('%0.3f',log10(et$table$PValue)),
 sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]
colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")

# Order by log fold change
o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]

# Save table
write.table(mat, file="DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")

#To exit R type the following
quit(save="no")
