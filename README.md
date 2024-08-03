
# RNA seq analysis of zika virus infected cell lines:

## Introduction

The Zika virus (ZIKV) has recently re-emerged in India, raising substantial concerns due to its association with microcephaly and other severe outcomes. Despite these alarming implications, the precise molecular mechanisms driving ZIKV pathogenesis are still not fully understood. This gap in knowledge underscores the need for further research to elucidate these mechanisms and identify potential biomarkers and genetic variants that could inform therapeutic strategies.

This project focuses on investigating ZIKV pathogenesis across different host cell lines. Our primary goal is to identify key biomarkers and decode the underlying genetic mechanisms contributing to ZIKV infection. Through an extensive analysis of RNA sequencing data, we have pinpointed several differentially expressed genes (DEGs) that reflect dysregulation in critical biological processes, including viral genome replication, the innate immune response, and metabolic pathways.

Our pathway analysis has highlighted the involvement of multiple signaling pathways—such as cytosolic DNA-sensing, Toll-like receptor, RIG-I-like receptor, TNF, and chemokine signaling pathways—in ZIKV pathogenesis. These findings suggest that these pathways play significant roles in the virus's interaction with host cells.

Additionally, the construction of a protein-protein interaction (PPI) network has identified ten hub genes (IFIT1, CXCL10, IFIT3, CXCL11, RSAD2, HERC5, ISG20, IFNB1, OAS1) as potential biomarkers for ZIKV infection. These genes are crucial for immune response, interferon regulation, viral replication, and neuronal damage, indicating their relevance to ZIKV-related diseases such as microcephaly.

The study also explored the impact of single nucleotide polymorphisms (SNPs) on infected cell lines. Notably, we identified potentially damaging SNPs within the GBP5 gene and other genes associated with ZIKV infection. These insights enhance our understanding of ZIKV pathogenesis and offer valuable information for early detection, risk assessment, and the development of targeted interventions.

This repository contains the data, methods, and results from our research, providing a comprehensive resource for those interested in ZIKV pathogenesis and related therapeutic strategies.

## Data Processing and Analysis Workflow

### 1. Downloading the Dataset from ENA Browser

**Purpose:** Obtain RNA-Seq data for analysis.

**Tool Used:** ENA Browser.

**Process:**
- Navigate to the [ENA Browser](https://www.ebi.ac.uk/ena/browser/home).
- Search for the dataset of interest using accession numbers (e.g., SRPxxxxxx, SRXxxxxxx, SRRxxxxxx).
- Download the data files (e.g., FASTQ files).

### 2. Pre-Alignment Quality Control

**Purpose:** Assess the quality of raw sequencing reads to ensure reliable data for subsequent analysis.

**Tool Used:** FastQC.

**Process:**
1. **Installation and Execution**:
   ```bash
   wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
   unzip fastqc_v0.11.9.zip
   cd FastQC
   chmod +x fastqc
   ```

2. **Quality Check**:
   ```bash
   cd data/raw
   fastqc *.fastq.gz
   xdg-open filename_fastqc.html
   ```
![image](https://github.com/user-attachments/assets/39073d05-1b9c-4674-8613-7586597b3fd5)

**Results:** FastQC provides metrics on read quality, including per-base sequence quality, GC content, and adapter content. Identify any issues such as low-quality reads or adapter contamination that need addressing before further analysis.

### 3. Trimming the Reads

**Purpose:** Remove adapters and low-quality bases to improve alignment accuracy.

**Tool Used:** Trimmomatic.

**Process:**
1. **Installation**:
   ```bash
   wget http://www.usadellab.org/cms/uploads/supported/Trimmomatic-0.39.tar.gz
   tar -xzvf Trimmomatic-0.39.tar.gz
   cd Trimmomatic-0.39
   ```

2. **Trimming**:
   ```bash
   java -jar trimmomatic-0.39.jar PE -phred33 input_R1.fastq.gz input_R2.fastq.gz trimmed_R1_paired.fastq.gz trimmed_R1_unpaired.fastq.gz trimmed_R2_paired.fastq.gz trimmed_R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
   ```

**Results:** Trimmed reads are cleaner, with removed adapters and improved quality. This step prepares the data for accurate alignment.

### 4. Alignment to Reference Genome

**Purpose:** Align RNA-Seq reads to the reference genome to identify the location of reads.

**Tool Used:** HISAT2.

**Process:**
1. **Installation**:
   ```bash
   wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2/download?path=%2F&files=hisat2-2.2.1-Linux_x86_64.zip
   unzip hisat2-2.2.1-Linux_x86_64.zip
   cd hisat2-2.2.1
   ```

2. **Prepare Reference**:
   ```bash
   wget ftp://ftp.ensembl.org/pub/release104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
   gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz
   hisat2-build Homo_sapiens.GRCh38.dna.toplevel.fa GRCh38_index
   ```

3. **Run Alignment**:
   ```bash
   hisat2 -x GRCh38_index -1 trimmed_R1_paired.fastq.gz -2 trimmed_R2_paired.fastq.gz -S aligned_reads.sam
   ```

4. **Convert and Sort BAM Files**:
   ```bash
   samtools view -bS aligned_reads.sam | samtools sort -o aligned_reads.sorted.bam
   samtools index aligned_reads.sorted.bam
   ```

**Results:** The alignment step produces BAM files with reads mapped to the reference genome, which are essential for downstream analyses.

### 5. Variant Calling

**Purpose:** Identify genetic variants from aligned RNA-Seq data.

**Tool Used:** BCFtools.

**Process:**
1. **Installation**:
   ```bash
   wget https://github.com/samtools/bcftools/releases/download/1.16.1/bcftools-1.16.1.tar.bz2
   tar -xjf bcftools-1.16.1.tar.bz2
   cd bcftools-1.16.1
   ```

2. **Call Variants**:
   ```bash
   bcftools mpileup -Ou -f reference_genome.fa aligned_reads.sorted.bam | bcftools call -mv -Ov -o variants.vcf
   ```

**Results:** The VCF file contains identified variants (SNPs and indels) from the RNA-Seq data, which can be further analyzed for their potential impact.

### 6. Annotating Variants with WANNOVAR

**Purpose:** Annotate variants to predict their functional impacts and associations.

**Tool Used:** WANNOVAR.

**Process:**
1. **Download and Installation**: Follow instructions from the [WANNOVAR website](http://wannovar.usc.edu/).

2. **Annotate Variants**:
   ```bash
   perl WANNOVAR.pl -i variants.vcf -o variants_annotated.vcf -d /path/to/annotation/files
   ```

**Results:** Annotated VCF file includes additional information such as predicted functional effects, gene associations, and potential clinical significance of each variant.

### 7. Visualizing Variants with VEP

**Purpose:** Further annotate and visualize the impact of variants using Ensembl's VEP tool.

**Tool Used:** VEP (Variant Effect Predictor).

**Process:**
1. **Install VEP**: Follow instructions from the [Ensembl VEP website](https://www.ensembl.org/info/docs/tools/vep/index.html).

2. **Run VEP**:
   ```bash
   vep -i variants_annotated.vcf --cache --dir_cache /path/to/cache --offline -o variants_vep.vcf
   ```

**Results:** VEP provides detailed annotations, including impact on genes and potential phenotypic effects. This enriches the understanding of variant consequences.

### 8. Gene Expression Analysis:


## Step 1: Run `featureCounts` - Quantification

In this step, you will use `featureCounts` to count reads mapped to genes from your BAM files. Follow the instructions below:

### 1.1 Download Annotation File

First, download the GTF annotation file if you haven't already:

```bash
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
```

### 1.2 Run `featureCounts`

Run `featureCounts` to quantify reads mapped to genes. Adjust the command below with the correct paths to your BAM files and the desired output directory:

```bash
featureCounts -p -a Homo_sapiens.GRCh38.106.gtf -o featurecounts/featurecounts6.txt aligned_read.sorted.bam
```

**Options:**
- `-p`: Counts only uniquely mapped reads.
- `-a`: Specifies the annotation file.
- `-o`: Specifies the output file.

### 1.3 Confirm Completion

After `featureCounts` finishes running, print a confirmation message:

```bash
echo "featureCounts finished running!"
```

### 1.4 Measure Execution Time (Optional)

To measure the execution time, you can use the following script:

```bash
SECONDS=0
# Run featureCounts command here
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
```
Note:convert the featurecounts.txt file to csv format  and prepare metadata file for differential expression analysis. 
## DESeq2 Analysis and Visualization

This script performs differential expression analysis using DESeq2 and visualizes the results with volcano and PCA plots.

### Prerequisites

Make sure you have the following R packages installed:
- `DESeq2`
- `ggplot2`

### Code

```r
library(DESeq2)
library(ggplot2)

# Read the count data and metadata from CSV files
count_data <- read.csv("C:\\Users\\ASUS\\Desktop\\annu\\assignment\\counts3.csv", header=TRUE, row.names=1)
meta_data <- read.csv("C:\\Users\\ASUS\\Desktop\\annu\\assignment\\meta3.csv", header=TRUE, row.names=1)

# Check if the sample names match between count_data and meta_data
if (!all(rownames(meta_data) %in% colnames(count_data))) {
  stop("Sample names in metadata and count data do not match.")
}

# Ensure the column names of count_data match the row names of meta_data
count_data <- count_data[, rownames(meta_data)]

# Check for NA values in the count_data
if (any(is.na(count_data))) {
  warning("NA values detected in count data. Replacing NA values with 0.")
  count_data[is.na(count_data)] <- 0
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta_data, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Save DESeq2 results to CSV
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# Extract upregulated and downregulated genes
upregulated_genes <- subset(res, log2FoldChange > 0 & padj < 0.05)
downregulated_genes <- subset(res, log2FoldChange < 0 & padj < 0.05)

# Save upregulated and downregulated genes to CSV files
write.csv(as.data.frame(upregulated_genes), file = "upregulated_genes.csv")
write.csv(as.data.frame(downregulated_genes), file = "downregulated_genes.csv")

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df$logP <- -log10(res_df$pvalue)
res_df$category <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 0, "Upregulated",
                          ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < 0, "Downregulated", "Not Significant"))

# Volcano plot with upregulated genes in blue and downregulated in red
ggplot(res_df, aes(x = log2FoldChange, y = logP, color = category)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not Significant" = "gray", "Upregulated" = "blue", "Downregulated" = "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

# Save the volcano plot
ggsave("volcano_plot.png", width = 8, height = 6)

# Perform variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Perform PCA
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate percent variance explained
percentVar <- attr(pca_data, "percentVar")

# PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(percentVar[1], 2), "% variance"),
       y = paste0("PC2: ", round(percentVar[2], 2), "% variance"),
       title = "PCA Plot") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

# Save the PCA plot
ggsave("pca_plot.png", width = 8, height = 6)

```

## 9. Protein-Protein Interaction Network Analysis

**Purpose:** Identify interactions between proteins encoded by the differentially expressed genes to understand their functional relationships.

**Tool Used:** STRING Database and Cytoscape.

**Process:**

1. **Access STRING:**
   - Go to [STRING Database](https://string-db.org/).
   - Input hub genes (e.g., IFIT1, CXCL10) to generate and visualize the PPI network.
   - Export the network data in a suitable format (e.g., TSV).

2. **Visualize with Cytoscape:**
   - Download and install Cytoscape from [Cytoscape.org](https://cytoscape.org/download.html).
   - Start Cytoscape and import the network data file (File > Import > Network > File).
   - Configure network visualization by adjusting node and edge styles, layouts, and labels.
   - Analyze network centrality metrics using NetworkAnalyzer (Tools > Network Analysis > NetworkAnalyzer) to identify hub genes.

**Results:** The STRING database provides a network of protein interactions, and Cytoscape enables detailed visualization and analysis. Identifying hub genes (highly connected nodes) can reveal key proteins involved in ZIKV infection.

## 10. Pathway Enrichment Analysis

**Purpose:** Identify biological pathways and processes significantly enriched with differentially expressed genes.

**Tool Used:** KEGG and DAVID.

**Process:**

1. **KEGG Analysis:**
   - Visit [KEGG Pathway](https://www.genome.jp/kegg/pathway.html).
   - Obtain KEGG Pathway IDs for your genes if needed.
   - Use the [KEGG Mapper](https://www.genome.jp/kegg/mapper.html) to search and visualize pathways associated with your gene list.

2. **DAVID Analysis:**
   - Go to [DAVID Bioinformatics Resources](https://david.ncifcrf.gov/).
   - Prepare a list of gene identifiers (e.g., gene symbols).
   - Upload the gene list (Upload > Paste or Upload).
   - Select `Functional Annotation Chart` and choose `KEGG Pathway` for pathway analysis.
   - Review and export the results in various formats (CSV, Excel).

**Results:** KEGG and DAVID analyses provide insights into the pathways significantly affected by ZIKV infection, helping to understand the biological context and functional implications of the differentially expressed genes.


## Conclusion

This repository contains detailed workflows and scripts for processing RNA-Seq data, performing differential expression analysis, variant calling, PPI network analysis with Cytoscape, and pathway enrichment using KEGG and DAVID. By following these steps, researchers can gain comprehensive insights into the molecular mechanisms of ZIKV infection and potentially identify new targets for therapeutic intervention.

   
