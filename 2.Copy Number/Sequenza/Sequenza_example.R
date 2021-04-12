# Install package
install.packages("sequenza")

# Load package
library(sequenza)
library(tidyverse)

# In the package is provided a small seqz file, we can use that (it is obtained
# from a BAM, see the Sequenza website if your are interested)
data.file <-
  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")

# The main interface consists of 3 functions:
#
# - sequenza.extract: process seqz data, normalization and segmentation
# - sequenza.fit: run grid-search approach to estimate cellularity and ploidy
# - sequenza.results: write files and plots using suggested or selected solution
#
# which run in sequence.

# Extract the information from the object
test <- sequenza.extract(data.file, verbose = FALSE)

# Example BAF plot
BAFs = Reduce(bind_rows,
              lapply(test$BAF %>% names, function(x)
                test$BAF[[x]] %>% mutate(chr = x))) %>% as_tibble()

BAFs %>% 
  ggplot(aes(x = start, y = mean)) +
  geom_point() +
  facet_wrap(~chr, scales = 'free_x')

# Sequenza fit (example chromosomes 1 and 10 used for the model, 1 is large)
CP <- sequenza.fit(test, chromosome.list = c(1,10))
CP

require(pheatmap)

pheatmap(
  CP$lpp,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color = RColorBrewer::brewer.pal(n = 8, "Reds")
)

# Write results
sequenza.results(sequenza.extract = test,
                 cp.table = CP, 
                 sample.id = "Test_Sample_Sequenza",
                 out.dir="Test_Sample_Sequenza")

# Loading mutation data
mutations = readr::read_tsv('Test_Sample_Sequenza/Test_Sample_Sequenza_mutations.txt')

ggplot(mutations,
       aes(`F`, fill = CNt %>% paste)) +
  geom_histogram(binwidth = 0.01)

# # More details on what happens inside Sequenza - GC normalisation (not shown in lecture
# which means you can skip this)
# 
# # GC content
# gc.stats <- gc.sample.stats(file = data.file, verbose = TRUE, parallel = 1)
# 
# gc_normal = gc.stats$normal
# gc_tumour = gc.stats$tumor
# 
# # Normalisation of the observed values by GC content
# normal = apply(
#   gc_normal$n,
#   1,
#   FUN = function(x, w) {
#     weighted.mean(x = w, w = x, na.rm = TRUE)
#   },
#   w = gc_normal$depth
# )
# 
# tumour = apply(
#   gc_tumour$n,
#   1,
#   FUN = function(x, w) {
#     weighted.mean(x = w, w = x, na.rm = TRUE)
#   },
#   w = gc_tumour$depth
# )
# 
# # Average tumour and normal depths
# avg_tum_depth <- weighted.mean(x = gc_tumour$depth, 
#                                w = colSums(gc_tumour$n))
# 
# avg_nor_depth <- weighted.mean(x = gc_normal$depth, 
#                                w = colSums(gc_normal$n))
