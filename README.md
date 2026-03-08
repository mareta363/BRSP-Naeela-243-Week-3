# Differential Expression Analysis of Immune Response to Influenza Vaccine

## Project Overview
This project investigates the **transcriptional immune response to influenza vaccination** using differential gene expression analysis of whole-blood transcriptomic data.

The goal of this analysis is to identify **genes that are differentially expressed after vaccination compared to baseline samples** and to interpret the biological pathways involved in the immune response.

The study uses publicly available microarray data and applies bioinformatics analysis in **R** to identify differentially expressed genes (DEGs) and perform functional enrichment analysis.

---

## Dataset

- **Dataset ID:** GSE48018  
- **Source:** NCBI Gene Expression Omnibus (GEO)  
- **Platform:** Illumina HumanHT-12 V4.0 Expression BeadChip (GPL10558)  
- **Samples:** Whole-blood RNA from healthy male volunteers  
- **Study:** Bucasas et al., 2011  

Baseline samples were collected **before vaccination**, and follow-up samples were collected **after influenza vaccination**.

GEO Links:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48018
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558

---

## Objectives

The objectives of this project are to:

- Identify **Differentially Expressed Genes (DEGs)** following influenza vaccination
- Visualize transcriptional changes using **volcano plots and heatmaps**
- Perform **functional enrichment analysis** to identify biological pathways
- Interpret immune mechanisms activated after vaccination

---

## Methods

### Data Retrieval
The dataset was retrieved from the **GEO database** using the `GEOquery` package in R.

### Data Preprocessing
The data underwent several preprocessing steps including:

- Quality control assessment
- **Log2 transformation**
- Inspection of normalization
- Evaluation of expression distribution using:
  - Boxplots
  - Histograms
  - UMAP visualization

### Differential Expression Analysis
Differential expression analysis was performed using the **limma** package based on linear modeling.

### Visualization
To interpret the results, several visualization methods were used:

- **Volcano plot** to identify significant DEGs
- **Heatmap** of the top 50 DEGs
- **UMAP** to examine sample clustering

### Functional Enrichment Analysis
Functional interpretation of DEGs was conducted using:

- **g:Profiler** for Gene Ontology (GO) enrichment analysis
- **KEGG Mapper** for pathway visualization

---

## Results

### Differential Gene Expression
The volcano plot revealed a **moderate number of upregulated genes and fewer downregulated genes** following vaccination. This pattern is consistent with vaccine-induced immune activation, which typically causes **targeted transcriptional changes rather than widespread gene expression alterations**.

Top upregulated genes included:

- **STAT1**
- **STAT2**
- **IFIT2**
- **GBP1**
- **OAS2**
- **TRIM22**

These genes are known **interferon-stimulated genes involved in antiviral immune responses**.

### Heatmap Visualization
A heatmap of the **top 50 DEGs** demonstrated gene expression differences between baseline and post-vaccination samples. Although global expression patterns remained similar, several immune-related genes showed **coordinated regulation**, indicating pathway-level immune activation.

---

## Functional Enrichment Analysis

### Gene Ontology (GO)

GO enrichment analysis identified several immune-related biological processes, including:

- **Defense response to another organism**
- **Antigen processing and presentation**
- **Cytolysis in another organism**

Cellular component enrichment included:

- **Cytosol**
- **ISGF3 complex** (STAT1–STAT2–IRF9 complex involved in interferon signaling)

Molecular function enrichment included:

- **Ribonucleotide binding**
- **Protein binding**

### KEGG Pathway Analysis

KEGG pathway mapping highlighted several important genes within the **Influenza A pathway**, including:

- **TLR7**
- **TRAF3**
- **TBK1**
- **STAT1**
- **STAT2**
- **OAS**
- **CASP1**

These genes participate in key immune mechanisms such as:

- Viral recognition
- Type I interferon signaling
- Activation of interferon-stimulated genes
- Inflammasome activation
- Antigen processing

---

## Key Findings

The analysis demonstrates that influenza vaccination induces a **coordinated antiviral immune response**, primarily mediated through:

- **Innate immune sensing**
- **Type I interferon signaling**
- **Activation of interferon-stimulated genes**

Although global gene expression patterns remained relatively stable, **specific immune pathways showed significant activation**, highlighting the targeted nature of vaccine-induced immune responses.

---

## Tools and Software

**Programming Environment**
- R
- RStudio

**R Packages**
- GEOquery
- limma
- ggplot2
- pheatmap
- umap

**Functional Analysis Tools**
- g:Profiler
- KEGG Mapper

---

## Author

**Naeela Mareta Ayudhia**  
BRSP Participant Number: **243**

---

## Reference

Bucasas, K. L., Franco, L. M., Shaw, C. A., Bray, M. S., Wells, J. M., Niño, D., Arden, N., Quarles, J. M., Couch, R. B., & Belmont, J. W. (2011).  
Early patterns of gene expression correlate with the humoral immune response to influenza vaccination in humans.  
*Journal of Infectious Diseases*, 203(7), 921–929.  
https://doi.org/10.1093/infdis/jiq156
