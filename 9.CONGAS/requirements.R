# for the heatmaps
install.packages("pheatmap")

# python interface
install.packages("reticulate")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# DE
BiocManager::install("DESeq2")

# gene annotations
BiocManager::install("biomaRt")

devtools::install_github("militeee/Rcongas")

# it may ask to install another instance of miniconda
# if you don't want problems just say yes (it just installa another
# instance of miniconda, it doe not remove other installations)
# in case you can also sincronize the working space
# with your current conda installation and environment

reticulate::py_install("congas==0.33", pip = T)