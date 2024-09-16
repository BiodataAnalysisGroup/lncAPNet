# lncAPNet

## Abstract

<div align='justify'> Transcriptomic analyses elucidate complexities in diseases like Chronic Lymphocytic Leukemia (CLL), revealing specific molecules, cellular traits, and signalling patterns influencing clinical outcomes. The status of the variable region of IGHV in leukemia cells is crucial for classifying CLL cases, providing insights into their clinical condition and prognosis. Patients can be stratified based on the mutational status of the IGHV genes, with mutated CLL (M-CLL) exhibiting a more indolent course than the more aggressive unmutated CLL (U-CLL) cases. Long non-coding RNAs (lncRNAs) play diverse roles in cancer but are often overlooked due to limitations in computational approaches. To address this, we propose lncAPNet [long non-coding Activity PASNet], a biologically informed Deep Learning model-based framework for supervised patient classification, instructed by lncRNA-mRNA-pathways interactions with eXplainable Artificial Intelligence (XAI) read-outs. This pipeline utilizes NetBID2 for network-based integration of omics data, and PASNet [Pathway Associated Sparse deep neural Network] for the deep learning analysis, further building on lncrnalyzer and EnrichR knowledge graphs as input for the biological priors. Complementing these, we deploy custom code integrating SHAP and heatmaps for interpretability, by combining a bespoke graph architecture representing lncRNA-mRNA-pathway interactions, facilitating a comprehensive analysis of regulatory mechanisms. The proposed framework is being tested by applying it into two extensive CLL transcriptomic datasets; data from ICGC-CLLES project (257 samples) are used as the training dataset, with the produced model validated on the BloodCancerMultiOmics2017 dataset (113 samples). Preliminary results showed that U-CLL cases are easily differentiated from M-CLL cases as expected, additionally revealing both biological ground-truth pathways, as well as insofar unknown but potential prognostic lncRNAs in disease. </div>

## Poster

![ECCB2024](https://raw.githubusercontent.com/BiodataAnalysisGroup/lncAPNet/main/Images/ECCB2024_poster.png)
