# Install package
install.packages("sequenza")

# Load package
library(sequenza)

# In the package is provided a small seqz file
data.file <-
  system.file("extdata", "example.seqz.txt.gz", package = "sequenza")

# The main interface consists of 3 functions:
#
# - sequenza.extract: process seqz data, normalization and segmentation
# - sequenza.fit: run grid-search approach to estimate cellularity and ploidy
# - sequenza.results: write files and plots using suggested or selected solution
#
# run them in sequence.

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

# Sequenza fit (example chromosome 10 used only for the model)
CP <- sequenza.fit(test, chromosome.list = 10)

# Write results
sequenza.results(sequenza.extract = test,
                 cp.table = CP, sample.id = "Test",
                 out.dir="TEST")

# More details on what happens inside Sequenza

# GC content
gc.stats <- gc.sample.stats(file = data.file, verbose = TRUE, parallel = 1)

gc_normal = gc.stats$normal
gc_tumour = gc.stats$tumor

# Normalisation of the observed values by GC content
normal = apply(
  gc_normal$n,
  1,
  FUN = function(x, w) {
    weighted.mean(x = w, w = x, na.rm = TRUE)
  },
  w = gc_normal$depth
)

tumour = apply(
  gc_tumour$n,
  1,
  FUN = function(x, w) {
    weighted.mean(x = w, w = x, na.rm = TRUE)
  },
  w = gc_tumour$depth
)

# Average tumour and normal depths
avg_tum_depth <- weighted.mean(x = gc_tumour$depth, 
                               w = colSums(gc_tumour$n))

avg_nor_depth <- weighted.mean(x = gc_normal$depth, 
                               w = colSums(gc_normal$n))
