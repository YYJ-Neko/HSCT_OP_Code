# HSCT_OP_Code
Analysis scripts for the hematopoietic stem cell transplantation associated microbiome study.

## Overview
This repository stores de-identified and streamlined code for statistical, microbiome and machine learning analyses.
All sensitive labels, local file paths and raw datasets have been removed. Raw clinical and sequencing data are not included.

The analysis pipeline covers:
- Microbial diversity analysis
- Differential biomarker screening
- Correlation and network analysis
- Metabolic function enrichment analysis
- Propensity score matching
- Survival and competing risk regression
- Regression model and visualization
- Python-based machine learning modeling

## Environment
R ≥ 4.4.1
Python ≥ 3.11.7

Main libraries:
cmprsk, survival, MatchIt, vegan, microeco, linkET, ReporterScore, igraph, tidyverse, ggplot2, scikit-learn, numpy, pandas

## Instructions
Prepare input data matching the structure required by each script.
Run R scripts for statistical and microbiome analyses.
Run Python scripts for machine learning modeling.
Visualization modules can be enabled on demand.

## Remarks
These scripts are provided for academic reproduction.
If part of the code is adopted, please cite the associated publication.
