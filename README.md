# BAGE
A Bayesian mixed linear model for age prediction based on PBMC gene expression data

## Overview
This repository contains all R code, workflows and resources associated with the research article:

> BAGE: A Bayesian framework for age prediction based on PBMC gene expression data
> 
> Veronica Suaste, Maria L. Daza-Torres, J. Cricelio Montesino-Lopez, Hilde Loge Nilsen.
> 
> Journal

In this research we used a Bayesian framework for chronological age prediction based on transcriptomic data (RNA-seq).
Particularly this analysis was developed with peripheral blood mononuclear cells(PBMCs) data.

## Data
All data used in the research is contained in the `./data` directory.

- __file1_healthy_samples.tsv: __

  Contains the metadata of all 174 healthy samples used in the research.
  The columns are metadata either found directly in the SRA, in the article related to the dataset or otherwise created from the original data.
  
  - id - identifier of the sequencing run data in the SRA (SRRXXX...)
  condition - either "healthy" or "disease" depending on whether the sample was control or patient respectively in the related study. In this case all samples have "healthy" condition, as this was the filtering factor.
  gender - either "male" or "female". This information was taken from SRA metadata or the related article supplementary files.
  - age - chronological age of the individual in years. This information was taken from SRA or the related article supplementary files.
  - library - either "paired" or "single"
  - instrument - sequencing instrument used to generate the reads. Information taken from metadata in the SRA.

- __counts.tsv :__ 

  Contains raw counts for the 174 samples that we used in the research. Each column is a sample named with the SRA identifier (SRRXXXX). 
  Raw counts were obtained using STAR aligner, see methods in article.
  
- __counts_rlog._tsv :__

  Contains the raw counts under the rlog transformation (see methods in article) for the 174 samples used in the research. Each column is a sample named with the SRA identifier (SRRXXXX). 

## Scripts

- run_models.R


## Reproduce article plots 

Data for generating plots is contained in the `./data/plots_data`

- __article_plots.Rmd :__  This notebook reproduce all plots in the article.
