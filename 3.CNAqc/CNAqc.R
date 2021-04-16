# Make it work in the folder where you downloaded the PCAWG data
# setwd('~/Downloads/WGS')

require(dplyr)

# PCAWG sample 
x = readRDS('0bfd1043-817c-e3e4-e050-11ac0c4860c5.rds')

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
# install.packages("devtools")
devtools::install_github("caravagnalab/CNAqc")

require(CNAqc)

cnaqc = CNAqc::init(
  snvs = x$mutations, # mutations
  cna = x$cna,        # copy number alterations
  purity = x$metadata$purity, # \pi (% of tumour cells in the biopsy)
  ref = 'hg19'        # Version of the reference genome we use (GRCh37)
)

cnaqc

# Plot segments
plot_segments(cnaqc)

# Plot data
plot_data_histogram(cnaqc, which = 'VAF')

# Peak detection to validate purity and ploidy inference
cnaqc = analyze_peaks(cnaqc, matching_epsilon = 0.025) 
cnaqc %>% plot_peaks_analysis(empty_plot = FALSE)

# Vignettes: https://caravagnalab.github.io/CNAqc/articles/a3_peaks_detection.html

# Results table
cnaqc$peaks_analysis$matches

# Overall score: this is a linear combination of the scores computed per expected-peak.
# The combination is computed for all the lines of the above table, weighting the 
# offset parameter by the weight coefficient. 
cnaqc$peaks_analysis$score

# The overall score is then used to determine the sample-level QC score (based on some
# maximum error tolerance matching_epsilon = 0.025 set above which is 5% overall tolerance)
cnaqc$peaks_analysis$QC

# Go and find as many PCAWG samples that have miscalled (wrongly called) values for
# ploidy/ purity of the input CNA segments!


