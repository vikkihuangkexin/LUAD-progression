# Inferring Evolutionary Trajectories from Cross-sectional Transcriptomic Data to Mirror Lung Adenocarcinoma Progression

# Introduction
Lung adenocarcinoma (LUAD) is a deadly tumor that remains approximately 15% survival rate in 5 years. The diverse reasons that made cancer progression and metastasis include genetic mutation, gene expression alteration, and so on. Understanding this dynamic process and identifying pivotal molecular events driving tumor progression is essential for improving LUAD diagnosis and treatment. Studying time-series data can simulate cancer evolution and determine the temporal patterns of molecular alterations. Unfortunately, it is difficult to collect complete time-series data from individual patients due to various reasons. In this study, we developed a trajectory model based on gene expression of LUAD patients and identified distinct progression branches. Moreover, we found a key gene BUB1B that could lead to aberrant spindle assembly checkpoint (SAC) signaling, and this appeared to be one of the key factors promoting LUAD progression. The trajectory shows clear evident mutation accumulation and clonal expansion along with the LUAD progression. We also investigated the clonal architectures and identified distinct clones and subclones in different LUAD branches. This clones and subclones can be used to identify combination therapy for LUAD. Our results demonstrate the biological utility and clinical application prospects of this progression model.


![Figure r21](https://user-images.githubusercontent.com/112677142/228104261-39291975-ab7f-4311-a833-3c06f84a372d.png)
 
A cartoon illustrating the reversed graph embedding method. Each point xi in the left panel refers a sample in the input space. yi refers this sample in the latent space. Each yi can map to the point fξ (yi ) by maintaining the graph structure through the reversed graph embedding. ε refers the connections between the points. By using this method, we can infer the underlying structure of high dimensional data.


<<<<<<< HEAD

Text below are brief introduction to the project code.

## Code Structure
This project is designed to provide a tool for lung cancer trajectory analysis. 

The main toy model code is located at `/LUAD/code/model/construction/LUAD_model.m`, which can be used to perform lung cancer trajectory analysis. The `/LUAD/code/model/construction/mRMR_0.9_compiled` and `/LUAD/code/model/construction/pami2016-pgsl-master` directories contain the compiled code for the feature selection function `mrmr_miq_d` and the original C++ code, respectively. Here we focus on explaining the `LUAD_model` code.

Dependencies
In order to run the code, you will need to have the following dependencies installed or put these code together with the model in the same folder:

* MATLAB
* mRMR 0.9 compiled code at `/LUAD/code/model/construction/mRMR_0.9_compiled`
* pami2016-pgsl-master code at `/LUAD/code/model/construction/pami2016-pgsl-master`

The toy model `LUAD_model.m` is well annotated, the model can be separated into these part:

* MRMR method feature selection
* PCA visualization
* Set parameters for the principal curve
* Calculate principal curve
* Drawing of principal curve
* Drawing of data
* Project data sample into principal curve

And the code in `/LUAD_model/code/function_analysis/` is the corresponding downstream analysis writen in R.

## Data needed for toy model

If you want to use our toy model easily, please modify the corresponding paths in the code to obtain the results in the paper. All data used in the code is available in the `/LUAD_model/data/` folder.


=======
# Usage
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9

