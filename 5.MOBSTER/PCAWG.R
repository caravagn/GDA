require(dplyr)

x = readRDS('~/Downloads/WGS/0168a2a6-c3af-4d58-a51c-d33f0fc7876d.rds')

cnaqc_x = CNAqc::init(
  x$mutations,
  x$cna,
  x$metadata$purity,
  'hg19'
)

# Diploid SNVs for mobster deconvolution
diploid_mutations = cnaqc_x %>%
  CNAqc::subset_by_segment_karyotype(karyotypes = '1:1') %>%
  CNAqc::subset_snvs()

diploid_mutations %>% CNAqc::plot_data_histogram()

x_fit = mobster::mobster_fit(diploid_mutations$snvs, auto_setup = 'FAST')
plot(x_fit$best)

