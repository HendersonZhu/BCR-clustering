# B-cell receptor (BCR) and T-cell receptor (TCR) repertoire clustering
*Method to cluster some-what biologically similar B-cell and T-cell receptors by Hamming distance in R.*

## *Introduction*
The vast B-cell repertoire data retrieved from RNA-seq would seem extremely diverse. It is convenient when biologically similar BCR can be classified and grouped as single entities to facilitate downsteam analysis. Hereby we proposed a method to cluster BCR with some-what similar CDR H3 amino acid sequences based on Hamming distance.

## *Requirements*
Load packages in R. Load any other metadata files required. Outputs from [VDJtools](https://github.com/mikessh/vdjtools) should be used as inputs for this clustering method (Shugay et al., 2015). \
R packages required for BCR clustering:
```
library("edgeR")
library("plyr")
library("dplyr")
library("data.table")
library("rlist")
library("textshape")
library("sjmisc")
library("reshape")
library("ggpubr")
library("Biostrings")
library("gdata")
library("stringr")
library("gridExtra")
library("tidyverse")
library("stringdist")
library("plotly")
library("seqinr")
```
To run clustering, please refer to /scripts/CDRH3_IGHV_convergence_clustering.R
## *Example of outputs*
Clusters associated to an IGHV gene between participants of different groups can be displayed:
![Image](/assets/CDRH3_cluster_expression_sample_1.png)
## *Project information*
This project is created and developed by Henderson Zhu (for BCR) and Caleb Batley (for TCR) under the supervision of Dr. Daniel O'Connor and Dr. Irina Chelysheva.
## *References*
Shugay M et al. VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. PLoS Comp Biol 2015; 11(11):e1004503-e1004503.
