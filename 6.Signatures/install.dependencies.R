# data.table
install.packages("data.table")

# BSgenome.Hsapiens.1000genomes.hs37d5
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# NMF
install.packages("NMF")

# nnls
install.packages("nnls")

# SparseSignatures
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SparseSignatures")

# lsa
install.packages("lsa")

# ggplot2
install.packages("ggplot2")

# gridExtra
install.packages("gridExtra")

# You need to install also SigProfilerExtractorR. The instructions to
# install it are available at https://github.com/AlexandrovLab/SigProfilerExtractorR. 
# 
# If your OS gives you trubles installing those, you can use conda.

