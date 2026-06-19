# lncAPNet: Deciphering lncRNA–mRNA Connections in Patient Transcriptomic Data

[![DOI](https://zenodo.org)](https://doi.org)

[![License: CC BY 4.0](https://shields.io)](https://creativecommons.org)

**lncAPNet** (long non-coding Activity PASNet) is an extended version of the APNet workflow designed for explainable, network-based patient stratification and supervised clustering using transcriptomic data. It integrates graph-based nonlinear inference of lncRNA–mRNA interactions with PASNet, a biologically informed sparse deep learning model.

---

## 📌 Motivation & Results

**lncAPNet** addresses the need for interpreting non-linear transcriptomic interactions in cancer, validated across:

* **CLL** (Chronic Lymphocytic Leukemia)

* **PRAD** (Prostate Adenocarcinoma)

* **BRCA** (Breast Invasive Carcinoma)

---

## 💾 Data Availability

Data for case studies is available on Zenodo [10.5281/zenodo.20637947](https://doi.org). Download the lncAPNet_v3.zip (2.7 GB) for necessary data matrices and pre-computed network files.

---

## 🛠️ Repository Structure

* /R/ — SJARACNe network generation and NetBID2/scMINER activity logic.

* /Python/ — PASNet deep learning architecture.

* /Nextflow/ — Automated execution pipelines.

---

## 🚀 Getting Started

1. **Clone:** git clone https://github.com

2. **Data:** Download/extract lncAPNet_v3.zip from Zenodo.

3. **Setup:** Install R (with NetBID2, scMINER), Python (PyTorch), and Nextflow.

---

## ✍️ Authors & License

* **Main Authors:** Vasileios Vasileiou, George Gavriilidis, Pedro Zeni, Marek Mraz, Evangelos Karatzas, Antonis Giakountis, Georgios Pavlopoulos, Antonis Giannakakis, Fotis Psomopoulos.

* **License:** [CC-BY-4.0](https://creativecommons.org)

---

## 🗺️ Citation

```bibtex
@dataset{vasileiou_2026_20637947,
  author       = {Vasileiou, Vasileios and Gavriilidis, George and Zeni, Pedro and Mraz, Marek and Karatzas, Evangelos and Giakountis, Antonis and Pavlopoulos, Georgios and Giannakakis, Antonis and Psomopoulos, Fotis},
  title        = {{lncAPNet enables the deciphering of lncRNA–mRNA connections in patient transcriptomic data}},
  publisher    = {Zenodo},
  year         = {2026},
  month        = {jun},
  version      = {v3},
  doi          = {10.5281/zenodo.20637947},
  url          = {https://doi.org}
}
```
