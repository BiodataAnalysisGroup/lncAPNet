\# lncAPNet: Deciphering lncRNA–mRNA Connections in Patient Transcriptomic Data



\[!\[DOI](https://zenodo.org)](https://doi.org)

\[!\[License: CC BY 4.0](https://shields.io)](https://creativecommons.org)



\*\*lncAPNet\*\* (long non-coding Activity PASNet) is an extended version of the APNet workflow designed for explainable, network-based patient stratification and supervised clustering using transcriptomic data. It integrates graph-based nonlinear inference of lncRNA–mRNA interactions with PASNet, a biologically informed sparse deep learning model.



\---



\## 📌 Motivation \& Results



\*\*lncAPNet\*\* addresses the need for interpreting non-linear transcriptomic interactions in cancer, validated across:

\*   \*\*CLL\*\* (Chronic Lymphocytic Leukemia)

\*   \*\*PRAD\*\* (Prostate Adenocarcinoma)

\*   \*\*BRCA\*\* (Breast Invasive Carcinoma)



\---



\## 💾 Data Availability



Data for case studies is available on Zenodo \[10.5281/zenodo.20637947](https://doi.org). Download the `lncAPNet\_v3.zip` (2.7 GB) for necessary data matrices and pre-computed network files.



\---



\## 🛠️ Repository Structure



\*   `/R/` — SJARACNe network generation and NetBID2/scMINER activity logic.

\*   `/Python/` — PASNet deep learning architecture.

\*   `/Nextflow/` — Automated execution pipelines.



\---



\## 🚀 Getting Started



1\.  \*\*Clone:\*\* `git clone https://github.com`

2\.  \*\*Data:\*\* Download/extract `lncAPNet\_v3.zip` from Zenodo.

3\.  \*\*Setup:\*\* Install R (with `NetBID2`, `scMINER`), Python (`PyTorch`), and Nextflow.



\---



\## ✍️ Authors \& License



\*   \*\*Main Authors:\*\* Vasileios Vasileiou, George Gavriilidis, Pedro Zeni, Marek Mraz, Evangelos Karatzas, Antonis Giakountis, Georgios Pavlopoulos, Antonis Giannakakis, Fotis Psomopoulos.

\*   \*\*License:\*\* \[CC-BY-4.0](https://creativecommons.org)



\---



\## 🗺️ Citation



```bibtex

@dataset{vasileiou\_2026\_20637947,

&#x20; author       = {Vasileiou, Vasileios et al.},

&#x20; title        = {{lncAPNet enables the deciphering of lncRNA–mRNA connections in patient transcriptomic data}},

&#x20; publisher    = {Zenodo},

&#x20; doi          = {10.5281/zenodo.20637947},

&#x20; url          = {https://doi.org}

}

```



