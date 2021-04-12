require(dplyr)
require(ggplot2)

####### Load all the input data, per sample, which is available in the
# samples folder. When you load it, add a column with the sample name to
# the loaded dataset

####### Load this table with the association between mutations and tumour type,
# then join the dataset loaded above with tumour type information

####### Filter out only clonal CNA segments using as reference the total_cn column,
# for instance just remove segments with > 6 copies

# Play around with this dataset and try to make the following visualisation
# - total number of segments per tumour type
# - per chromosome barplot of total CN per tumour type
# - barplot as above, but retaining only diploid segments (1 = Major, 1 = minor)
# - overall (whole-genome) plot of 
# - rank chromosomes by number or size of LOH events (0 = minor) 


# Extra: make a genome-wide segmentation plot like the one we have seen for
# Sequenza, drawing minor and major segments across the whole genome. In this
# case we need to transform relative coordinates into absolute coordinates.
# Absolute ones depend on the overall size of chromosomes, which I provide you with.
# 
# First, install this package
# devtools::install_github('caravagnalab/CNAqc')
# 
# Second, get this tibble with coordinates
# CNAqc::chr_coordinates_hg19
# 
# Write your adjustment and visualisation function.


  
  
  