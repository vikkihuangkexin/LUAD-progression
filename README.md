# Inferring Evolutionary Trajectories from Cross-sectional Transcriptomic Data to Mirror Lung Adenocarcinoma Progression

# Introduction
Lung adenocarcinoma (LUAD) is a deadly tumor that remains approximately 15% survival rate in 5 years. The diverse reasons that made cancer progression and metastasis include genetic mutation, gene expression alteration, and so on. Understanding this dynamic process and identifying pivotal molecular events driving tumor progression is essential for improving LUAD diagnosis and treatment. Studying time-series data can simulate cancer evolution and determine the temporal patterns of molecular alterations. Unfortunately, it is difficult to collect complete time-series data from individual patients due to various reasons. In this study, we developed a trajectory model based on gene expression of LUAD patients and identified distinct progression branches. Moreover, we found a key gene BUB1B that could lead to aberrant spindle assembly checkpoint (SAC) signaling, and this appeared to be one of the key factors promoting LUAD progression. The trajectory shows clear evident mutation accumulation and clonal expansion along with the LUAD progression. We also investigated the clonal architectures and identified distinct clones and subclones in different LUAD branches. This clones and subclones can be used to identify combination therapy for LUAD. Our results demonstrate the biological utility and clinical application prospects of this progression model.


![Figure r21](https://user-images.githubusercontent.com/112677142/228104261-39291975-ab7f-4311-a833-3c06f84a372d.png)
 
A cartoon illustrating the reversed graph embedding method. Each point xi in the left panel refers a sample in the input space. yi refers this sample in the latent space. Each yi can map to the point fξ (yi ) by maintaining the graph structure through the reversed graph embedding. ε refers the connections between the points. By using this method, we can infer the underlying structure of high dimensional data.


# Usage

