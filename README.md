# iJRF

## Introduction
The integrative Joint Random Forest (iJRF) algorithm, developed by Petralia et al. (Bioinformatics, 2017), is a machine-learning framework designed to infer gene regulatory networks across multiple related experimental conditions simultaneously. Unlike traditional random forest–based approaches that model each condition independently, iJRF constructs joint forests in which the same regulatory features (e.g., miRNAs) are evaluated across all conditions when splitting nodes, thereby promoting the identification of both shared and condition-specific regulatory interactions. The method models the expression of each mRNA as a function of multiple candidate regulators using random forests and introduces an integrative mechanism that borrows statistical strength across conditions, improving robustness in small-sample, high-dimensional settings. iJRF further incorporates biological prior knowledge—such as known miRNA–mRNA relationships—by biasing the selection of predictors during tree construction, ensuring that the inferred networks are both data-driven and biologically informed. The output consists of a set of condition-specific networks that capture common and differential regulatory structures, enabling the discovery of molecular rewiring events associated with environmental exposures, disease states, or other perturbations.

## Citation
For more information, please visit or cite these articles: 

Petralia, Francesca, Pei Wang, Jialiang Yang, and Zhidong Tu. "Integrative random forest for gene regulatory network inference." Bioinformatics 31, no. 12 (2015): i197-i205.

Petralia, Francesca, Won-Min Song, Zhidong Tu, and Pei Wang. "New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer." Journal of proteome research 15, no. 3 (2016): 743-754.

Petralia, F., Aushev, V.N., Gopalakrishnan, K., Kappil, M., W Khin, N., Chen, J., Teitelbaum, S.L. and Wang, P., 2017. A new method to study the change of miRNA–mRNA interactions due to environmental exposures. Bioinformatics, 33(14), pp.i199-i207.


## Installing iJRF R package

* library(devtools)
* install_github("petraf01/iJRF")
* install_github("WangLab-MSSM/BayesDeBulk/BayesDeBulk")
