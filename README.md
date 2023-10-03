# DLBCL_CoExprLandscape

This repository includes the code necessary to reproduce the results from the paper "The co-expression landscape of diffuse large B-cell lymphoma". For the ARACNe
implementation used for coexpression values calculation please refer to https://github.com/ddiannae/ARACNE-multicore

- 1_DataProcessing: This scripts downloads TCGA RNA-seq data from DLBCL samples, performs pre-processing, DEG analysis and normalization of the data.
- 2_FunctionalEnrichment : This script performs functional enrichment analysis on a list of communities (see Supplementary material for formatting examples). 
- 3_BPSummary : Integrates the functional enrichment of a network with a set of DEG analysis.
- 4_SurivalAnalysis: Performs Kaplan-Meier analysis on the DLBCL samples using the median expression value for a set of given genes.
