# Make it work in the folder where you downloaded the PCAWG data
# setwd('~/Downloads/WGS')

require(dplyr)

# PCAWG sample 
x = readRDS('00db1b95-8ca3-4cc4-bb46-6b8c8019a7c7.rds')

x %>% names

# Inspect what is in these fields
x$mutations %>% 
  select(chr, from, to, ref, alt, DP, NV, VAF, Hugo_Symbol, Variant_Classification, everything()) %>% 
  View

# CNAs
x$cna %>% 
  View

# Metadata 
x$metadata

# We now work together with
# 
# CNAqc - https://caravagnalab.github.io/CNAqc/
# 
# devtools::install_github("caravagnalab/CNAqc")
