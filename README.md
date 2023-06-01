
# AlphaDiv

<!-- badges: start -->
<!-- badges: end -->

This repo provides the code and example data to reproduce the results in the manuscript **"A Mixed-Effect Kernel Machine Regression Model for Integrative Analysis of Alpha Diversity in Microbiome Studies"**.

## Overview
Our method is motivated by the need to integrate data from multiple microbiome studies for reproducible discovery of the association between microbiome alpha diversity and a phenotype of intereset. Individual studies often provide inconsistent results due to insufficient sample size, heterogeneous study populations and study variability. 

In this paper, we propose a mixed-effect kernel machine regression model to assess the association of alpha diversity with a phenotype of interest. Our approach readily incorporates the study-specific characteristics (including sequencing protocols) to allow for flexible modeling of the microbiome effect via a kernel similarity matrix. Within the proposed framework, we provide three hypothesis testing approaches to answer research questions of interest, including (1) common effect (2) heterogeneity-by-protocol (HBP) effect (3) overall effect.

## Instructions
To run the code, you need to clone the repo.
```
git clone https://github.com/stephlee3/AlphaDiv.git
cd AlphaDiv
```

Install and load all the necessary packages.

```{r}
library(tidyverse)
library(nlme)
library(lme4)
library(expm)
library(Matrix)
library(CompQuadForm)
```

Now you can run the example code in `example.R`. P-values from (1) common effect (2) heterogeneity-by-protocol (HBP) effect (3) overall effect will be generated on the simulated data. 

## Basic Structure
This repo is organized as follows

* `Code/`
  * `common_effect.R`: the function for the test of common effect.
  * `HBP_effect.R`: the function for the test of HBP effect.
  * `overall_effect.R`: the function for the test of overall effect.
  
* `Simulation/`
  * `model_based.R`: the reproducible code for model-based simulation.
  * `community_based.R`: the reproducbile code for community-based simulation.

* `example.R`: the example code to run the tests on simulated data.


## Contact Information
For any additional questions, please contact Runzhe Li(rli51@jhmi.edu).




