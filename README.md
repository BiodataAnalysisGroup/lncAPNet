# lncAPNet

## Abstract

<div align='justify'> Transcriptomic analyses elucidate complexities in diseases like Chronic Lymphocytic Leukemia (CLL), revealing specific molecules, cellular traits, and signalling patterns influencing clinical outcomes. The status of the variable region of IGHV in leukemia cells is crucial for classifying CLL cases, providing insights into their clinical condition and prognosis. Patients can be stratified based on the mutational status of the IGHV genes, with mutated CLL (M-CLL) exhibiting a more indolent course than the more aggressive unmutated CLL (U-CLL) cases. Long non-coding RNAs (lncRNAs) play diverse roles in cancer but are often overlooked due to limitations in computational approaches. To address this, we propose lncAPNet [long non-coding Activity PASNet], a biologically informed Deep Learning model-based framework for supervised patient classification, instructed by lncRNA-mRNA-pathways interactions with eXplainable Artificial Intelligence (XAI) read-outs. This pipeline utilizes NetBID2 for network-based integration of omics data, and PASNet [Pathway Associated Sparse deep neural Network] for the deep learning analysis, further building on lncrnalyzer and EnrichR knowledge graphs as input for the biological priors. Complementing these, we deploy custom code integrating SHAP for interpretability, by combining a bespoke graph architecture representing lncRNA-mRNA-pathway interactions, facilitating a comprehensive analysis of regulatory mechanisms. The proposed framework is being tested by applying it into two extensive CLL transcriptomic datasets; data from ICGC-CLLES project (257 samples) are used as the training-validation dataset, with the produced model tested on the BloodCancerMultiOmics2017 dataset (113 samples). Preliminary results showed that U-CLL cases are easily differentiated from M-CLL cases as expected, additionally revealing both biological ground-truth pathways, as well as insofar unknown but potential prognostic lncRNAs in disease. </div>

## Poster

![ECCB2024](https://raw.githubusercontent.com/BiodataAnalysisGroup/lncAPNet/main/Images/ECCB2024_poster.png)

## Framework steps

**lncAPNet framework on the [ICGC](https://dcc.icgc.org/) and [BCMO](https://www.bioconductor.org/packages/release/data/experiment/html/BloodCancerMultiOmics2017.html) datasets aims to identify lncRNAs that act as molecular drivers**

* **Step1:**
<br><br>
**Run [SJARACNe](https://github.com/jyyulab/SJARACNe)/[NetBID2](https://jyyulab.github.io/NetBID/) for Network Reconstruction and Differential Activity Analysis**

**a.** Prepare and preprocess the ICGC and BCMO datasets, focusing on both **mRNA** and **lncRNA** expression data.

**b.** Use **SJARACNe** to reconstruct **gene regulatory networks (GRNs)**, identifying interactions between mRNAs and lncRNAs.

**c.** Run **NetBID2** to estimate the activity levels of genes and lncRNAs as **driver molecules**.

**d.** Perform **Differential Activity analysis** to highlight lncRNAs and mRNAs that are potential drivers influencing key biological pathways, especially focusing on lncRNAs as potential master regulators.

* **Step2:**
<br><br>
**Match Pathways with Drivers and Use PASNet for Interpretability**
  
**a.** Match identified drivers (including lncRNAs and mRNAs) to biological pathways using pathway databases [EnrichR-KG](https://maayanlab.cloud/enrichr-kg), [lncRNAlyzr](https://github.com/MaayanLab/lncRNAlyzr) (***Gene Ontology, KEGG, Reactome, Wikipathways***).

**b.** Run **[PASNet](https://github.com/DataX-JieHao/PASNet)** to integrate biological pathway information and train a deep learning model on the ICGC dataset (***70-30 split for training/validation***), and then *test* on the BCMO dataset.

**c.** Utilize **[SHAP values](https://github.com/shap/shap)** for deep learning explainability, identifying which lncRNAs and pathways are most important in the model’s predictions.


* **Step3:**
<br><br>
**Post-hoc Network Analysis for Biological Explainability**
  
**a.** Perform a **Post-hoc network analysis** to explore the connections between identified lncRNA/mRNA drivers and their related pathways.

**b.** Use network visualization tools to map these relationships, focusing on the role of lncRNAs as drivers.

**c.** Interpret the network to uncover how lncRNAs function as **key regulators** of biological processes.
