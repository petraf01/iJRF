# iJRF

* [Introduction](#introduction)
* [Citation](#citation)
* [Installing iJRF package in R]


## Introduction
The integrative Joint Random Forest (iJRF) algorithm, developed by Petralia et al. (PMID: 26072483, PMID: 26072483 ), provides an integrated framework that includes implementations for both iRafNet and JRF.

* IRAFNet (PMID: 26072483) is a computational framework designed to integrate multi-omic data—such as transcriptomic, proteomic, and knock-out experiments to infer gene regulatory networks. By combining probabilistic modeling with prior biological knowledge, IRAFNet identifies key regulatory modules and pathways driving specific phenotypes or disease states.

* JRF (Joint Random Forest) is a statistical framework developed to jointly infer co-expression networks across multiple related biological conditions or data sets. Unlike traditional single-condition models, JRF leverages shared information across contexts to improve network reconstruction accuracy while allowing for condition-specific regulatory differences. It employs a joint modeling approach that encourages both sparsity and similarity across inferred networks, capturing conserved and differential gene–gene interactions. This enables the identification of associations that are stable across conditions as well as those that uniquely characterize specific disease states, tissue types, or experimental perturbations.

## Citation
For more information, please visit or cite these articles: 

* Petralia, Francesca, Pei Wang, Jialiang Yang, and Zhidong Tu. "Integrative random forest for gene regulatory network inference." Bioinformatics 31, no. 12 (2015): i197-i205.

* Petralia, Francesca, Won-Min Song, Zhidong Tu, and Pei Wang. "New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer." Journal of proteome research 15, no. 3 (2016): 743-754.

* Petralia, F., Aushev, V.N., Gopalakrishnan, K., Kappil, M., W Khin, N., Chen, J., Teitelbaum, S.L. and Wang, P., 2017. A new method to study the change of miRNA–mRNA interactions due to environmental exposures. Bioinformatics, 33(14), pp.i199-i207.


## Installing iJRF R package

* library(devtools)
* install_github("petraf01/iJRF")
