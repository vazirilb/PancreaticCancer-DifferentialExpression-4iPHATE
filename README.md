# PancreaticCancer-DifferentialExpression-4iPHATE

Repository for HPNE-KRAS(G12V) RNA-seq, TCGA-PAAD analysis and CellPhate assessment using iterative indirect immunofluorescence imaging (4i) performed for manuscript titled 'TRIP13 protects pancreatic cancer cells against intrinsic and therapy-induced DNA replication stress' submitted to NAR Cancer.

## Assessing KRAS-induced DNA damage resposne (DDR) transcriptome:

### RNA-sequencing of hTERT-HPNE cells expressing KRASG12V
#### RNA isolation, library preparation and sequencing
Total RNA was isolated from hTERT-immortalized HPNE-DT cells with or without doxycycline induction of KRASG12V for six hours using QIAGEN RNeasy (QIAGEN, 74104). Total RNA was quantified using Invitrogen Qubit, ran on an agarose gel to check RNA integrity and were submitted to Novogene (Sacramento, CA) for library preparation and sequencing. Novogene performed preliminary quality control (QC), sample quantification and sample integrity test using agarose gel electrophoresis, Nanodrop and Agilent 2100 bioanalyzer. All the samples passed the quality control tests. Library preparation and sequencing were performed by Novogene. cDNA libraries containing 250-300 bp cDNA inserts were enriched for mRNA using poly(A)-selection and sequenced using the Illumina PE150 platform.
#### RNA-seq data preprocessing
RNAseq analysis pipeline is available in the ‘HPNE RNASeq analysis.R’ file. Briefly, raw paired-end sequencing reads in the form of FASTQ files were assessed for quality control and trimmed using Trim Galore. Trimmed reads were aligned to the GRCh38.p12 human reference genome using STAR. Transcript-level abundances were estimated using Salmon in alignment-based mode. Transcript-level estimates were summarized to the gene level using the tximport package. Count matrices were rounded and loaded into a DGEList object using edgeR. Genes were filtered to retain only those classified as protein-coding and located on autosomal chromosomes (excluding those on sex/mitochondrial chromosomes). Filtered gene-level count matrices were exported as CSV files for downstream differential expression analysis and visualization.
#### RNA-seq data analysis
DESeq2 (v1.42.1) was used to generate differentially expressed genes (DEGs) for pathway analysis in R (v4.3.3). Gene Set Enrichment Analysis (GSEA) was performed using the `clusterProfiler` (v4.10.1) package to identify enriched Gene Ontology (GO) terms for biological processes (BP) in the differentially expressed genes. A variance stabilizing transformation (VST) performed for principal component analysis (PCA) (plotPCA function) and hetamaps (pheatmap, v1.0.12). 

### TCGA-PAAD analysis
The TCGA pancreatic adenocarcinoma (TCGA-PAAD) data was retrieved using TCGAbiolinks R package (v2.30.0). Samples defined as “high purity” and “treatment-naïve” were used for gene expression analysis. Data analysis and visualization were performed as described above for hTERT-HPNE samples. Survival analysis and Kaplan-Meier curves were generated using the original TCGA_PAAD data and was performed using the following R packages: survival (v3.2-13) and survminer (v0.4.9).


