RNA-Sequence-Analysis
A comprehensive RNA-Seq analysis pipeline designed to process and analyze transcriptomic data from different conditions (e.g., disease vs. normal or UHR vs. HBR samples). This README outlines the workflow from raw data acquisition to differential expression analysis, mutation annotation, and gene set enrichment analysis.

Table of Contents
Loading and Overview of Data
Quality Control
Adapter Trim
Alignment
SAM to BAM Conversion
Merge BAM Files
Expression
HTSEQ-Count
Differential Expression (Ballgown)
Differential Expression with HTSEQ-Count (edgeR)
Visualization of DE Genes with Ballgown
Analysis and Visualization of DE Genes without Ballgown
RNA Sequence Mutation (ANNOVAR)
Gene Set Enrichment Analysis (GSEA)
Scripts Used
References and Links
Loading and Overview of Data
Commands to download raw data from SRA and convert .sra files to FASTQ format:

bash
Copy
# Download raw data from SRA
prefetch SRR28709957 SRR28709958 SRR28709959 SRR28709960 SRR28709961 SRR28709962 \
        SRR28709963 SRR28709964 SRR28709965 SRR28709966 SRR28709967 SRR28709968

# Convert downloaded .sra files to FASTQ format
fastq-dump --split-files SRR28709957.sra
fastq-dump --split-files SRR28709958.sra
# ... repeat for all SRR accession numbers ...
Quality Control
Commands used for quality control:

bash
Copy
# Perform quality control checks on raw FASTQ files
fastqc *.fastq

# Summarize FastQC reports across all samples using MultiQC
multiqc .
Adapter Trim
Commands to trim adapters using AfterQC (example for each SRR):

bash
Copy
python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709957_1.fastq.gz \
  -2 raw_data/SRR28709957_2.fastq.gz \
  > afterqc_SRR28709957.log 2>&1

python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709958_1.fastq.gz \
  -2 raw_data/SRR28709958_2.fastq.gz \
  > afterqc_SRR28709958.log 2>&1

python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709959_1.fastq.gz \
  -2 raw_data/SRR28709959_2.fastq.gz \
  > afterqc_SRR28709959.log 2>&1

# ... repeat for all SRR accession numbers ...
Upon completion of trimming, the cleaned sequences are stored in the good directory, e.g.:

python-repl
Copy
SRR28709957_1.good.fq.gz
SRR28709957_2.good.fq.gz
...
SRR28709968_1.good.fq.gz
SRR28709968_2.good.fq.gz
A total of 24 cleaned FASTQ files were generated (12 samples, paired-end reads).

Alignment
Commands to align using HISAT2 (example for each SRR):

bash
Copy
hisat2 -p 8 \
  --rg-id=SRR28709957 --rg SM:UHR --rg LB:SRR28709957 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 \
  --dta --rna-strandness RF \
  -x /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome/ \
  -1 /media/data07/sude/workspace/rnaseq/data/myProject/good/SRR28709957_1.good.fq.gz \
  -2 /media/data07/sude/workspace/rnaseq/data/myProject/good/SRR28709957_2.good.fq.gz \
  -S /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good.sam \
  --summary-file /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957.hisat2.output \
  > /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_hisat2.log 2>&1

# ... repeat for SRR28709958, SRR28709959, etc. ...
SAM to BAM Conversion
Commands to convert and sort SAM files to BAM using samtools:

bash
Copy
samtools sort -@ 8 \
  -o /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good.sam && \
samtools index -@ 8 /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam

samtools sort -@ 8 \
  -o /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good.sam && \
samtools index -@ 8 /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam

# ... repeat for all SRR accession numbers ...
Merge BAM Files
Use Picard’s MergeSamFiles to combine replicates (example):

bash
Copy
java -Xmx2g -jar /media/data07/sude/workspace/rnaseq/student_tools/picard.jar MergeSamFiles \
  OUTPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/UHR.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709959_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709960_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709961_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709962_good.bam

java -Xmx2g -jar /media/data07/sude/workspace/rnaseq/student_tools/picard.jar MergeSamFiles \
  OUTPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/HBR.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709963_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709964_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709965_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709966_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709967_good.bam \
  INPUT=/media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709968_good.bam
Expression
Download and prepare annotation (GTF):

bash
Copy
wget -O /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf.gz \
  ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

gunzip /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf.gz
StringTie transcript abundance estimation (example commands):

bash
Copy
stringtie -p 8 -G /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf \
  -e -B \
  -o UHR_Rep57/transcripts.gtf \
  -A UHR_Rep57/gene_abundances.tsv \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam

stringtie -p 8 -G /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf \
  -e -B \
  -o UHR_Rep58/transcripts.gtf \
  -A UHR_Rep58/gene_abundances.tsv \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam

# ... repeat for each sample ...
Merge all samples into a single expression matrix:

bash
Copy
./stringtie_expression_matrix.pl --expression_metric=TPM \
  --result_dirs='HBR_Rep63,HBR_Rep64,HBR_Rep65,HBR_Rep66,HBR_Rep67,HBR_Rep68,
                 UHR_Rep57,UHR_Rep58,UHR_Rep59,UHR_Rep60,UHR_Rep61,UHR_Rep62' \
  --transcript_matrix_file=transcript_tpm_all_samples.tsv \
  --gene_matrix_file=gene_tpm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=FPKM \
  --result_dirs='HBR_Rep63,HBR_Rep64,HBR_Rep65,HBR_Rep66,HBR_Rep67,HBR_Rep68,
                 UHR_Rep57,UHR_Rep58,UHR_Rep59,UHR_Rep60,UHR_Rep61,UHR_Rep62' \
  --transcript_matrix_file=transcript_fpkm_all_samples.tsv \
  --gene_matrix_file=gene_fpkm_all_samples.tsv

./stringtie_expression_matrix.pl --expression_metric=Coverage \
  --result_dirs='HBR_Rep63,HBR_Rep64,HBR_Rep65,HBR_Rep66,HBR_Rep67,HBR_Rep68,
                 UHR_Rep57,UHR_Rep58,UHR_Rep59,UHR_Rep60,UHR_Rep61,UHR_Rep62' \
  --transcript_matrix_file=transcript_coverage_all_samples.tsv \
  --gene_matrix_file=gene_coverage_all_samples.tsv
HTSEQ-Count
Commands to generate gene counts for each sample:

bash
Copy
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 \
  --type exon --idattr gene_id \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf \
  > UHR_Rep57_gene.tsv

htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 \
  --type exon --idattr gene_id \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf \
  > UHR_Rep58_gene.tsv

# ... repeat for all SRR accession numbers ...
Merge all HTSeq count files into a single table:

bash
Copy
join HBR_Rep63_gene.tsv | join - HBR_Rep64_gene.tsv | join - HBR_Rep65_gene.tsv | join - HBR_Rep66_gene.tsv \
| join - HBR_Rep67_gene.tsv | join - HBR_Rep68_gene.tsv | join - UHR_Rep57_gene.tsv | join - UHR_Rep58_gene.tsv \
| join - UHR_Rep59_gene.tsv | join - UHR_Rep60_gene.tsv | join - UHR_Rep61_gene.tsv | join - UHR_Rep62_gene.tsv \
> gene_counts_all.tsv

echo "GeneID HBR_Rep63 HBR_Rep64 HBR_Rep65 HBR_Rep66 HBR_Rep67 HBR_Rep68 \
      UHR_Rep57 UHR_Rep58 UHR_Rep59 UHR_Rep60 UHR_Rep61 UHR_Rep62" \
> header.txt

cat header.txt gene_counts_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' \
> gene_counts_final.tsv

rm -f gene_counts_all.tsv header.txt
Differential Expression (Ballgown)
Differential gene expression analysis is a critical step in RNA-Seq studies, aiming to identify genes with statistically significant expression changes between conditions. We utilized Ballgown to analyze normalized FPKM/TPM values for the UHR vs. HBR comparison.

Total genes analyzed: 62,710
Genes passing expression filter (variance > 1, etc.): 14,269
Differentially expressed genes (p-value < 0.05): 855
Below are the top DEGs with higher abundance in UHR and HBR. For detailed R code, see R_Part1Ballgown.R (as mentioned in English).

Differential Expression with HTSEQ-Count (edgeR)
For count-based RNA-Seq data, we used the edgeR Bioconductor package. Raw counts from HTSeq were input into edgeR to identify DEGs between UHR and HBR.

Total genes analyzed: 62,710
Genes passing expression filter: 14,269
DEGs (p-value < 0.05, FDR < 0.01): 855
We selected top DEGs based on log-fold change (log2FC) and statistical significance. Below is a Venn diagram showing the overlap of DEGs identified by Ballgown (left circle) and edgeR (right circle).

(Replace this link with the actual image path in your repo, e.g., images/Ballgown_edgeR_venn.png.)

Visualization of DE Genes with Ballgown
Here, we visualize and analyze RASD2, a gene known for its involvement in neurodevelopmental disorders and various cancers. RASD2 appeared in both Ballgown and HTSeq analyses, showing significantly high TPM in the intersecting DEGs between UHR and HBR.

Biological Significance of RASD2 in Cancer
RASD2 is crucial in signal transduction (GTPase activity), neurogenesis, and synaptic plasticity.
Linked to tumor progression and drug resistance via PI3K/AKT, ERK, and Ras pathways.
Overexpression correlates with poor prognosis and metastatic potential in certain cancers.
Identified in glioblastomas, lung cancer, and colorectal cancer, making it a potential biomarker.
Our analysis shows RASD2 is differentially expressed in UHR vs. HBR, with higher TPM values. See R_Part2Ballgown.R for code and detailed plots.

Analysis and Visualization of DE Genes without Ballgown
We also explored transcriptomic differences using FPKM values, correlation analysis, multidimensional scaling (MDS), and heatmap visualization (see Supplementary.R). This confirms significant differences between UHR and HBR, highlighting the biological relevance of our dataset. All images can be inspected for deeper insights.

RNA Sequence Mutation (ANNOVAR)
We performed variant calling and used ANNOVAR to annotate genetic variants from RNA-Seq data:

Variant calling using an RNA-seq mutation pipeline (rnaseqmut).
ANNOVAR mapped variants to genes and predicted their effects.
Initially, 76,922 variants were detected.
We filtered out non-coding, synonymous, and high-frequency variants in 1000 Genomes.
We prioritized exonic, nonsynonymous mutations and used multiple tools (PolyPhen-2, SIFT, MutationTaster, PROVEAN, CADD) to evaluate pathogenicity.
Finally, 20 significant mutations remained (clinically or functionally relevant).
More details and raw CSV are available here:
Google Drive Link

Gene Set Enrichment Analysis (GSEA)
GSEA determines whether predefined sets of genes show significant differences between two states. We performed GSEA on:

HTSeq-Ballgown Intersection Genes (top 10 DEGs).
ANNOVAR Variant Genes (top 20 mutated genes).
Pathways Analyzed
Gene Ontology (GO): Biological processes, molecular functions, cellular components.
KEGG: Metabolic and signaling pathways.
Immune Pathway: Immune-related gene interactions.
Key Findings
GO enrichment for the top 10 DEGs revealed significant immune and DNA metabolism processes (e.g., antigen receptor signaling, B cell receptor signaling, complement activation).
Immune response emerges as a key regulatory mechanism.
For details, see GeneSetEnrinchmentAnalysis.R (or similarly named script). Result CSV files and plots are here:
Google Drive Link

Scripts Used
R_Part1Ballgown.R
Ballgown DE analysis (filtering low-abundance genes, statistical testing).
R_Part2Ballgown.R
Ballgown-based visualization, gene-level plots (e.g., RASD2).
edgeR.R
Differential expression analysis with edgeR on HTSeq count data.
Supplementary.R
Additional visualizations (heatmaps, MDS, correlation plots) without Ballgown.
GeneSetEnrinchmentAnalysis.R
GSEA for top DEGs and mutated genes.
References and Links
Ballgown: https://www.bioconductor.org/packages/release/bioc/html/ballgown.html
edgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html
HTSeq: https://htseq.readthedocs.io/en/master/
StringTie: https://ccb.jhu.edu/software/stringtie/
ANNOVAR: http://annovar.openbioinformatics.org/en/latest/
GSEA: http://www.broadinstitute.org/gsea/
Google Drive (DEGs, Mutations, Enrichment Results):
https://drive.google.com/drive/folders/1Ty43qhz2SRW3w50ksGvpMQ8GD9ju3FR2?usp=sharing
End of README

Feel free to explore each step in detail, modify parameters as needed, and consult the scripts for deeper insight into the analysis workflow. If you have questions or suggestions, please open an issue or submit a pull request.
