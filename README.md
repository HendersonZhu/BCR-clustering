# B-cell receptor (BCR) repertoire clustering
*Method to cluster some-what biologically similar B-cell receptors by Hamming distance.*

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
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("Biostrings")
library("gdata")
library("stringr")
library("gridExtra")
library("tidyverse")
library("stringdist")
library("plotly")
library("seqinr")
```
## *Example of outputs*

## *References*
Shugay M et al. VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. PLoS Comp Biol 2015; 11(11):e1004503-e1004503.
