# lncAPNet

**lncAPNet** is an extension of the Activity PASNet (**APNet**) framework that incorporates **long non-coding RNAs (lncRNAs)** into explainable, network-based patient stratification and supervised clustering using transcriptomic data.

This repository hosts the implementation and resources for the **lncAPNet** computational pipeline.

---

## Overview

lncAPNet is a novel computational pipeline built upon the recently reported **Activity PASNet (APNet)** framework  
([Gavriilidis et al., 2025](https://doi.org/10.1093/bioinformatics/btaf063)).  
It performs **explainable AI (XAI)-driven supervised clustering of patient cohorts** by integrating:

- Transcriptomic analysis  
- Biological priors  
- The concept of **differential activity**

The original APNet framework integrates:
- **SJARACNe** for gene regulatory network inference  
- **NetBID2** for driver activity analysis  
- **PASNet**, a sparse and interpretable deep learning architecture for patient classification  

These components are supported by biological priors and post hoc visualization modules for clinical bioinformatics applications.

---

## Extensions Introduced by lncAPNet

lncAPNet expands the APNet toolbox by explicitly modeling **lncRNA-driven regulatory mechanisms**. The key extensions include:

- Graph-based, nonlinear interpretation of **lncRNA–mRNA regulatory relationships** using the SJARACNe / NetBID2 framework  
- Integration of an **lncRNA-specific knowledge graph (lncRNAlyzr-KG)** into the biological priors  
  (Evangelista et al., 2025)  
- Identification and interpretation of **lncRNA-associated regulons** derived from NetBID2 analysis  
- Enhanced explainability of lncRNA contributions to patient stratification and disease mechanisms  

---

## Pipeline Status

- ([Nextflow Installation instructions]([https://doi.org/10.1093/bioinformatics/btaf063](https://nf-co.re/docs/get_started/environment_setup/nextflow)))

- A **Nextflow-based implementation** of the lncAPNet pipeline is now available for **bulk RNAseq** to improve scalability, portability, and reproducibility.


🚧 **For scRNAseq data is under active development**

---

## Implementation

## **bioRxiv** link as pre-print: 
([Vasileiou V. et al., 2026](https://doi.org/10.64898/2025.12.18.695074))
