README
======

This repository contains three scripts to reproduce some of the results
published in the article: 

**Pascual-García, A., & Bell, T. Community-level signatures of ecological succession in natural bacterial communities. _Nat Commun_ 11, 2386 (2020). [DOI: 10.1038/s41467-020-16011-3](https://doi.org/10.1038/s41467-020-16011-3)**
  

## INPUT DATA ##

Input data can be found in [Zenodo](https://doi.org/10.5281/zenodo.3539537). Please dowload the input files and locate them in the folder `source`.

## SCRIPTS ##

There are three scripts available. Each script has a short description in the header and a number of options to fix between the statements _START EDIT_ and _STOP EDIT_. The scripts are ready to be executed. It is required the same structure of directories found in the repository, in which case only editing the root pathway for the repository is needed.

* **LocationRandHierarchy2method\_V2.R**

	Script to compute either the ANOSIM, ADONIS2 or MRPP metric for the observed data across different classifications in the spatial distance (Figure 1 in the MS)

* **LocationRandHierarchy2syntheticData.R**

   Script to generate synthetic data to understand the behaviour of the above metrics.
   
* **SEM_2pub.R**
   
   Script to compute structural equation models for the functional measurements. A number of models are provided in the folder `LavaanModels`. The rational of the organization of these folders is the following. In the first exploratory stage, we investigated the basic causal relationships between the variables and the effect of working with latent variables (folder _Exploratory Models_), leading to a first satisfactory structural model (Milestone1). In the second stage (folder _Milestone1_), we introduced the definition of the different classes which, given the large number of new parameters introduced, led to a worst model (as determined by the AIC value). Therefore, we reespecified the model systematically looking for modifications of this basic model (changes in the type of relationship or addition of new relations, labelled “mod”) or deletions of links (labelled “del”), leading to a second model (_Milestone2_). The reespecified model had few modifications, in particular the deletion of a pathway between P and Cells and the inclusion of all the correlations between exoenzymatic variables. The third stage investigates constraints between the variables in the partial linear coefficients (folders _Milestone2_ and _Milestone 3_). The models presented in one of the Suppl. Tables are available in the main folder and labelled as "free" (of constrains) "fixed" (i.e. all the parameters constrained to be the same across groups) and the final selected model which has no additional label. To generate the fourth model presented in the table simply select the final model and activate the option to shuffle the samples in the code.
   
It is also included a modified version of the R function heatmap.2 (gplots) that allows for changing the size of the labels and that is sourced in one of the scripts.

For more details, please refer to the original article. 

