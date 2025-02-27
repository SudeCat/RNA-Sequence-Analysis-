# RNA-Sequence-Analysis

A comprehensive RNA-Seq analysis pipeline designed to process and analyze transcriptomic data from different conditions (e.g., disease vs. normal or UHR vs. HBR samples). This repository outlines the workflow from raw data acquisition to differential expression analysis, mutation annotation, and gene set enrichment analysis.

- [Loading and Overview of Data](#loading-and-overview-of-data)
- [Quality Control](#quality-control)
- [Adapter Trim](#adapter-trim)
- [Alignment](#alignment)
- [SAM to BAM Conversion](#sam-to-bam-conversion)
- [Merge BAM Files](#merge-bam-files)
- [Expression](#expression)
- [HTSEQ-Count](#htseq-count)
- [Differential Expression (Ballgown)](#differential-expression-ballgown)
- [Differential Expression with HTSEQ-Count (edgeR)](#differential-expression-with-htseq-count-edger)
- [Visualization of DE Genes with Ballgown](#visualization-of-de-genes-with-ballgown)
- [Analysis and Visualization of DE Genes without Ballgown](#analysis-and-visualization-of-de-genes-without-ballgown)
- [RNA Sequence Mutation (ANNOVAR)](#rna-sequence-mutation-annovar)
- [Gene Set Enrichment Analysis (GSEA)](#gene-set-enrichment-analysis-gsea)



## Loading and Overview of Data
Commands to download raw data from SRA and convert .sra files to FASTQ format:

```bash
# Download raw data from SRA
prefetch SRR28709957 SRR28709958 SRR28709959 SRR28709960 SRR28709961 SRR28709962 \
        SRR28709963 SRR28709964 SRR28709965 SRR28709966 SRR28709967 SRR28709968

# Convert downloaded .sra files to FASTQ format
fastq-dump --split-files SRR28709957.sra
fastq-dump --split-files SRR28709958.sra
... repeat for all SRR accession numbers ... 
```

## Quality Control
Commands used for quality control:

```bash
# Perform quality control checks on raw FASTQ files
fastqc *.fastq

# Summarize FastQC reports across all samples using MultiQC
multiqc .
```
## Adapter Trim
Commands to trim adapters using AfterQC (example for each SRR):

```bash

python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709957_1.fastq.gz \
  -2 raw_data/SRR28709957_2.fastq.gz > afterqc_SRR28709957.log 2>&1

python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709958_1.fastq.gz \
  -2 raw_data/SRR28709958_2.fastq.gz > afterqc_SRR28709958.log 2>&1

python2.7 /media/data07/sude/workspace/rnaseq/student_tools/AfterQC/after.py \
  -1 raw_data/SRR28709959_1.fastq.gz \
  -2 raw_data/SRR28709959_2.fastq.gz > afterqc_SRR28709959.log 2>&1

# ... repeat for all SRR accession numbers ...
```
Upon completion of the trimming process, the cleaned sequences were stored in the good directory with the naming convention:

```bash
SRR28709957_1.good.fq.gz
SRR28709957_2.good.fq.gz
SRR28709958_1.good.fq.gz
SRR28709958_2.good.fq.gz
...
SRR28709968_1.good.fq.gz
SRR28709968_2.good.fq.gz
```
A total of 24 cleaned FASTQ files were generated (12 samples, paired-end reads).

## Alignment
Commands to align using HISAT2 (example for each SRR):

```bash

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
```
# SAM to BAM Conversion
Commands to convert and sort SAM files to BAM using samtools:

```bash
samtools sort -@ 8 -o /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good.sam && \
samtools index -@ 8 /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709957_good_sorted.bam

samtools sort -@ 8 -o /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam \
  /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good.sam && \
samtools index -@ 8 /media/data07/sude/workspace/rnaseq/data/myProject/alignment/SRR28709958_good_sorted.bam

# ... repeat for all SRR accession numbers ...
```

## Merge BAM Files
Picard’s MergeSamFiles can be used to combine replicates:

```bash
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
```
## Expression
Download and prepare annotation (GTF):

```bash
wget -O /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf.gz \
  ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

gunzip /media/data07/sude/workspace/rnaseq/data/myProject/alignment/grch38/genome.gtf.gz
```
StringTie transcript abundance estimation (example commands):

```bash
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
```
Merge all samples into a single expression matrix:

```bash
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
```
## HTSEQ-Count
Commands to generate gene counts for each sample:

```bash
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
```
Merge all HTSeq count files into a single table:

```bash
join HBR_Rep63_gene.tsv | join - HBR_Rep64_gene.tsv | join - HBR_Rep65_gene.tsv | join - HBR_Rep66_gene.tsv \
| join - HBR_Rep67_gene.tsv | join - HBR_Rep68_gene.tsv | join - UHR_Rep57_gene.tsv | join - UHR_Rep58_gene.tsv \
| join - UHR_Rep59_gene.tsv | join - UHR_Rep60_gene.tsv | join - UHR_Rep61_gene.tsv | join - UHR_Rep62_gene.tsv \
> gene_counts_all.tsv

echo "GeneID HBR_Rep63 HBR_Rep64 HBR_Rep65 HBR_Rep66 HBR_Rep67 HBR_Rep68 UHR_Rep57 UHR_Rep58 UHR_Rep59 UHR_Rep60 UHR_Rep61 UHR_Rep62" \
> header.txt

cat header.txt gene_counts_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_counts_final.tsv

rm -f gene_counts_all.tsv header.txt

```
## Differential Expression (Ballgown)
Differential gene expression analysis is a critical step in RNA-Seq studies, aiming to identify genes with statistically significant expression changes between conditions. In this project, I utilized Ballgown to perform differential expression analysis on normalized FPKM/TPM values.

Total genes analyzed: 62,710
Genes passing the expression filter (variance > 1, removal of low-abundance): 14,269
Differentially expressed genes (p-value < 0.05): 855
I focused on genes showing significant changes between UHR and HBR conditions. Below, the Top 20 DEGs with higher abundance in each condition were identified. For the detailed R code, see [R_Part1Ballgown.R](R_Part1Ballgown.R).


### Differential Expression with HTSEQ-Count (edgeR)
For count-based RNA-Seq data, I used the edgeR Bioconductor package. Raw counts from HTSeq were input into edgeR to identify differentially expressed genes (DEGs) between UHR and HBR. You can check the script here [edgeR.R](edgeR.R).

Total genes analyzed: 62,710
Genes passing expression filter: 14,269
DEGs (p-value < 0.05, FDR < 0.01): 855
I also selected top differentially expressed genes based on log-fold change (log2FC) and statistical significance. Genes with higher abundance in UHR and HBR were identified separately.

### Venn Diagram of Overlaps
Below is a Venn diagram showing the overlap of DEGs identified by Ballgown (left circle) and edgeR (right circle).

![venny](images/venny.png)


## Visualization of DE Genes with Ballgown
In this section, I visualize and analyze the expression levels of RASD2, a gene known for its involvement in neurodevelopmental disorders and various cancers. RASD2 appeared in both Ballgown and HTSeq differential expression analysis outputs, showing significantly high TPM values in intersecting DEGs between UHR and HBR conditions.

### Biological Significance of RASD2 in Cancer
RASD2 (Ras-related protein RASD2) is crucial in signal transduction pathways, especially those related to GTPase activity, neurogenesis, and synaptic plasticity.
It has been linked to tumor progression and drug resistance, interacting with oncogenic pathways like PI3K/AKT, ERK, and Ras signaling.
Overexpression in certain cancers correlates with poor prognosis and increased metastatic potential.
RASD2 has been identified in glioblastomas, lung cancer, and colorectal cancer, making it a potential biomarker for tumor progression.
My analysis suggests that RASD2 is differentially expressed in UHR vs. HBR conditions, with a notable increase in TPM values. See [R_Part2Ballgown.R](R_Part2Ballgown.R) for the code used for visualization and further analysis.

### Analysis and Visualization of DE Genes without Ballgown
Beyond Ballgown, I also explored transcriptomic differences using FPKM expression values, correlation analysis, multidimensional scaling (MDS), and heatmap visualization. Differentially expressed genes were highlighted and exported, confirming significant differences between UHR and HBR. All relevant R code is found in [Supplementary.R](Supplementary.R)
and the [images](images/)
 can be inspected for further insights.

## RNA Sequence Mutation (ANNOVAR)
Identifying genetic variants in RNA-Seq data is crucial for understanding disease pathogenesis and functional genomics. I used ANNOVAR to annotate variants from RNA sequencing data.

Variant calling on RNA-Seq data was performed via the rnaseqmut pipeline.
ANNOVAR mapped the identified variants to genes and predicted effects on protein function.
A total of 76,922 genetic variants were initially detected.
Filtering criteria included:
Retaining exonic and nonsynonymous mutations (removing synonymous and non-coding variants).
Excluding variants found at high frequency in 1000 Genomes Project.
Using multiple computational tools (PolyPhen-2, SIFT, MutationTaster, PROVEAN, CADD) to assess functional impact.
Retaining variants in ClinVar (disease-associated) and COSMIC (cancer-related).
After these filters, 20 significant mutations remained. These are likely pathogenic or disease-related. You can [check the raw CSV file and further details here](https://drive.google.com/drive/folders/1Ty43qhz2SRW3w50ksGvpMQ8GD9ju3FR2?usp=sharing ). 


## Gene Set Enrichment Analysis (GSEA)
Gene Set Enrichment Analysis (GSEA) identifies whether predefined sets of genes show significant, concordant differences between two states. I performed GSEA on:

### HTSeq-Ballgown Intersection Genes: The top 10 DEGs common to both Ballgown and edgeR.
#### ANNOVAR Variant Genes
The top 20 mutated genes identified from ANNOVAR analysis.

**Pathways Analyzed
Gene Ontology (GO)**: Biological processes, molecular functions, and cellular components.

**KEGG**: Known metabolic and signaling pathways.

**Immune Pathway**: Immune-related gene interactions.

**Key Findings**

GO Enrichment for the top 10 DEGs showed significant processes in immune response and DNA metabolism (e.g., antigen receptor signaling, B cell receptor signaling, humoral immune response).
The immune response appears to be a key regulatory mechanism in this dataset.
All relevant code for this section is in [GeneSetEnrichmentAnalysis.R](GeneSetEnrichmentAnalysis.R).
Result CSV files and additional plots can [be found here](https://drive.google.com/drive/folders/1Ty43qhz2SRW3w50ksGvpMQ8GD9ju3FR2).



## End of the README

Feel free to explore each step in detail, modify parameters as needed, and consult the scripts for deeper insight into the analysis workflow. If you have questions or suggestions, please open an issue or submit a pull request.
- **Email:** [catsudeebrar@gmail.com](mail)
