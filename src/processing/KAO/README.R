## README ##

## This folder contains processing files for the metabolomics data.

## "Time_stamp.R"
# description: This file writes the time stamp for metabolomics raw files;
#  the time stamp is used for run-order correction.
# date created: 4/5/2019
# date last modified: 4/5/2019

## "GC_runOrderCorrection_R_20200810.R"
# description: This file takes time stamp information and metabolite quant
#  values to correct for run-time effect (i.e. deviations in signal which 
#  modulate over time). The strategy is also to adjust the mean signal to 
#  be the same across batches by using the pooled wild type means ( see line 163).
#  Metabolite features are filtered for data quality (see line 213) based on 
#  feature 'tier'. Tier 1 is the best, Tier 5 indicates imputated quantification.
#  This file also generates exploratory data analysis plots (regressions and PCA) 
#  and outputs normalized data files, which were used for later analysis. 
# date created: 8/10/2020
# date last modified: 12/16/2020

## 