# Multi-Omics Analysis of Polar Organic Compounds in Daphnia magna

## 🔬 Aims of the Study

This MSc Bioinformatics group project investigates the impact of polar organic compounds extracted from natural river environments on *Daphnia magna* using advanced machine learning and integrative multi-omics approaches. The goals are to:

- Explore how modern machine learning (e.g. IRF, MOFA) can be applied to complex, real-world multi-omics data.
- Assess biological effects of polar organic compounds via in vivo bioassays.
- Integrate transcriptomic and metabolomic data to reveal hidden patterns and candidate biomarkers.

---

## 🧪 Experimental Design

### 1. River Organic Extracts
- Samples were collected from **12 sites (D01–D12)** along a river, enriched 100,000–250,000× over environmental levels.
- Extracts were split for:
  - Chemical analysis (91 semi-quantified compounds)
  - Bioassay exposure

### 2. Treatment Media
- Organic extracts were diluted to **1× and 10×** environmental levels in standard culture media with 0.08% methanol.
- Controls were standard media with methanol only.

### 3. 48-Hour Exposure Bioassay
- **15 neonates** of *Daphnia magna* per treatment, with **6 biological replicates** each.
- Exposed for 48 hours, no feeding.
- Samples frozen post-exposure for omics extraction.

### 4. Transcriptomics
- **150 samples** sequenced using BGI DNBseq platform.
- Data normalized via **DESeq2**.

### 5. Metabolomics
- Measured using **LTQ-Orbitrap Elite** in both positive and negative ion modes.
- Preprocessed with **DIMSpy**, PQN normalization, and glog transformation.

### 6. Feature Annotation
- Genes mapped to:
  - *Drosophila melanogaster*
  - *Homo sapiens*
- Metabolite peaks annotated using **BEAMSpy** and **KEGG/HMDB**:
  - KEGG
  - HMDB

---

## 🎯 Research Objectives

1. **How does omics data change** across locations (D01–D12), treatment concentrations (Control, 1×, 10×), and their interactions?
2. **What are the cross-omics interactions** and how do they change under these conditions?
3. **Which individual chemicals** are most associated with adverse outcomes in *Daphnia*?
4. What **biological or environmental insights** can be drawn from the patterns and features discovered?

---

## 🔍 Approach & Analysis Overview

To address the research objectives, the following analysis workflow was implemented:

### 1. Data Preprocessing
- **Outlier and missing sample removal**: Omics datasets were filtered to exclude samples with missing data or technical outliers, ensuring high-quality input for downstream analysis.
- **Chemical data augmentation**: Synthetic replicates were generated for chemical measurements to enable compatibility with statistical learning models and increase robustness.

### 2. Machine Learning: Single-Omics & Multi-Omics Analysis
- **Single-omics analysis**: Iterative Random Forests (IRF) were applied separately to transcriptomics and metabolomics datasets to:
  - Identify features (genes or metabolites) most predictive of site and treatment concentration.
  - Visualize classification performance using ROC curves.
  - Extract Gini-based importance scores for each feature.
  
- **Multi-omics integration**: Multi-Omics Factor Analysis (MOFA2) was used to:
  - Uncover latent factors driving variation across omics layers.
  - Visualize clustering and separation of samples by site and treatment.
  - Identify top contributing features to each latent factor (e.g., Factor 3).

### 3. Biological Interpretation
- **Feature annotation and ortholog mapping**:
  - Genes were mapped to human and *Drosophila* orthologs using OrthoDB.
  - Metabolite peaks were annotated using KEGG and HMDB databases.

- **Pathway enrichment**:
  - Top-ranked genes and metabolites (based on IRF and MOFA weights) were input to **HMDB and KEGG pathway databases** to identify enriched biological pathways.
  - These were interpreted in the context of the chemical exposure profiles across different river sites, revealing site-specific chemical signatures and associated biological processes.

This integrative analysis pipeline allowed us to link environmental chemical signatures with transcriptomic and metabolomic responses in *Daphnia magna*, and provided insights into potentially adverse outcome pathways at different geographic locations and treatment concentrations.

---

## 🧰 Computational Environment

This project was executed on a **High-Performance Computing (HPC)** environment using the following R versions:

- **R version 4.4.1**: Used for all steps including MOFA, data preprocessing, visualisation, and feature annotation.
- **R version 4.3.1**: Required specifically for the **Iterative Random Forest (IRF)** analysis due to package dependencies.

All R scripts were run in isolated environments using module loading to ensure reproducibility.

---

## 📁 Repository Structure

```text
project-root/
│
├── scripts/
│ ├── omics_data_preprocessing.R # Omics loading and preprocessing
│ ├── chemical_data_preprocessing.R # Omics loading and preprocessing
│ ├── IRF.R # Iterative Random Forest on omics (R 4.3.1)
│ ├── MOFA.R # MOFA model training and visualization
│ ├── ID_Mapping_RNA.R # Gene ID conversion
│ └── ID_Mapping_Metabolites.R # Metabolite ID conversion
│
├── figures/
│ ├── irf/ # IRF plots
│ ├── mofa/ # MOFA plots
│ └── preprocessing/ # preprocessing plots
│
└── output/
  └── mofa_model.hdf5 # Saved MOFA model
```