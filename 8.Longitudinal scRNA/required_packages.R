if (!require("progress")) install.packages("progress")
if (!require("parallel")) install.packages("parallel")
if (!require("Rfast")) install.packages("Rfast")
if (!require("Seurat")) install.packages('Seurat')
if (!require("dplyr")) install.packages('dplyr')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

if (!require("devtools")) install.packages("devtools")
devtools::install_github("BIMIB-DISCo/LACE", ref = "master")