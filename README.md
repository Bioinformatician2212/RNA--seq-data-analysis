
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

### 8. Gene Expression Analysis

**Purpose:** Identify differentially expressed genes between conditions to understand ZIKV infection effects.

**Tool Used:** DESeq2.

**Process:**
1. **Install DESeq2**:
   ```R
   install.packages("BiocManager")
   BiocManager::install("DESeq2")
   ```

2. **Generate Read Counts**:
   ```bash
   featureCounts -a annotation.gtf -o gene_counts.txt aligned_reads.sorted.bam
   ```

3. **Run DESeq2**:
   ```R
   library(DESeq2)
   count_data <- read.table("gene_counts.txt", header=TRUE, row.names=1)
   col_data <- data.frame(condition=factor(c("control", "treatment", "control")))
   dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
   dds <- DESeq(dds)
   res <- results(dds)
   write.csv(as.data.frame(res), file = "deseq2_results.csv")
   ```

4. **Create Volcano and PCA Plots**:
   ```R
   library(ggplot2)
   res <- read.csv("deseq2_results.csv", row.names=1)
   ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point(aes(color=padj < 0.05), alpha=0.5) + scale_color_manual(values=c("black", "red")) + theme_minimal() + labs(x="Log2 Fold Change", y="-Log10 p-value", color="Significant")
   ggsave("volcano_plot.png")

   rld <- rlog(dds, blind=FALSE)
   pcaData <- plotPCA(rld, intgroup = "condition", returnData=TRUE)
   ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size=3) + theme_minimal() + labs(x="PC1", y="PC2", color="Condition")
   ggsave("pca_plot.png")
   ```

5. **Extract Upregulated and Downregulated Genes**:
   ```R
   res <- read.csv("deseq2_results.csv", row.names=1)
   alpha <- 0.05
   log2fc_threshold <- 1
   upregulated <- subset(res, padj < alpha & log2FoldChange > log2fc_threshold)
   downregulated <- subset(res, padj < alpha & log2# RNA--seq-data-analysis

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

   
