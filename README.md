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
- Integration of an **lncRNA-specific knowledge graph (lncRNAlyzr-KG)** into the biological priors ([Evangelista et al., 2025](https://doi.org/10.1016/j.jmb.2025.168938))
- Identification and interpretation of **lncRNA-associated regulons** derived from NetBID2 analysis  
- Enhanced explainability of lncRNA contributions to patient stratification and disease mechanisms  

---

## Pipeline Status

### Nextflow - Main Pipeline

- [Nextflow Installation Instructions](https://nf-co.re/docs/get_started/environment_setup/nextflow)

- A **Nextflow-based implementation** of the lncAPNet pipeline is now available for **bulk RNAseq** to improve scalability, portability, and reproducibility.

### Post-Hoc Analysis Requirments

#### Requirments

Nextflow/Setup/PostHocAnalysis_requirments.txt

🚧 **Post-hoc analyses Python Colab notebooks (connected with Nextflow pipeline) is under active development**

- ✅ 1.PASNet_output_analysis.ipynb
  
- ✅ 2.PostHoc_Graphs.ipynb
  
- 🚧 3.Clinical_Covariates.ipynb
  
- ✅ 4.DataBases_Validation.ipynb

### scRNAseq data

🚧 **For scRNAseq data is under active development**

---

## Implementation

Define the local path for `input` in the `Nextflow/pipeline/conf/local.config` file.
If using singularity, in the same file, uncomment and update the local path to where the singularity containers are stored.

The pipeline should then run using `docker` or `singularity` containers with this command
from within the `Nextflow/pipeline` folder:
```
nextflow run main.nf -c ../conf/local.config -profile singularity,local --outdir output
```
or, directly passing the input to the CLI command:
```
nextflow run main.nf -profile docker --input /path/to/data/eset.rds --gene_info /path/to/data/gene_info.xlsx --group0 "GROUP0" --group1 "GROUP1" --comparison "COMPARISON" --outdir output
```

If you are using a Unix terminal, you need to convert the scripts to Unix format using:
```
dos2unix ./bin/*
```

### Example from Nextflow_demo/data toy dataset

```
nextflow run main.nf -profile docker --input ../Nextflow_demo/data/test_eset.rds --outdir output --gene_info ../Nextflow_demo/data/gene_info.xlsx --group0 "M" --group1 "U" --comparison "IGHV"
```

## **bioRxiv** link as pre-print: 
([Vasileiou V. et al., 2026](https://doi.org/10.64898/2025.12.18.695074))
