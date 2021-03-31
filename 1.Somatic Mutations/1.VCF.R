# Packages
require(dplyr)
require(ggplot2)

# install.packages("vcfR")
require(vcfR) 

# VCF - Variant Calling Format - of a multi-region tumour (WGS ~80x median coverage)
VCF_url = "https://raw.githubusercontent.com/caravagnalab/CNAqc_datasets/main/MSeq_Set06/Mutations/Set.06.WGS.merged_filtered.vcf"

# Download, load and cancel data
download.file(VCF_url, "Set.06.WGS.merged_filtered.vcf",)

set6 = vcfR::read.vcfR("Set.06.WGS.merged_filtered.vcf")

file.remove("Set.06.WGS.merged_filtered.vcf")

# VCF
print(set6)

# We extract all the information we need, using the tidy data representation format.

# INFO fields 
info_tidy = vcfR::extract_info_tidy(set6)

# Fixed fields (mutation coordinates, chr|from|ref|alt)
fix_tidy = set6@fix %>%
  as_tibble %>%
  rename(
    chr = CHROM,
    from = POS,
    ref = REF,
    alt = ALT
  ) %>%
  mutate(from = as.numeric(from), to = from + nchar(alt))

# Genotypes
geno_tidy = vcfR::extract_gt_tidy(set6) %>%
  group_split(Indiv)

# Sample mutations in the CNAqc format
sample_mutations = lapply(
  geno_tidy, 
  function(x) 
  {
    bind_cols(info_tidy, fix_tidy) %>%
      full_join(x, by = "Key") %>%
      mutate(DP = as.numeric(gt_NR), NV = as.numeric(gt_NV)) %>%
      mutate(VAF = NV / DP) %>%
      select(chr, from, to, ref, alt, NV, DP, VAF, everything()) %>%
      filter(!is.na(VAF), VAF > 0) # VAF > 0 in each sample
  })

# A list for all samples available
names(sample_mutations) = sapply(sample_mutations, function(x) x$Indiv[1])
sample_mutations = sample_mutations[!is.na(names(sample_mutations))]

# Design of these calls (from ICR/UCL, UK)
# 
# - multiple tumour samples for the same patient
# - one germline
# 
# 1) Use Mutect2 to call tumour vs normal, for each sample
# 2) Gather all called mutations (across all tumour-normal pairs)
# 3) Genotype alltoghether the mutations with Platypus (germline caller)
# 
# This helps calling low-frequency variants in one sample that have large
# support in other samples (they would be removed by a single-sample analysis,
# while they are called in a genotyping run).
# 
# We have output FILTERS from Platypus.
# 
# "alleleBias" makes no sense for tumour, it is for germline. It means the VAF is 
# not close than 0.5, which of course in tumour is due to aneuploidy and
# the intrinsic fact that somatic SNVs are not germline SNPs. So the FILTER flags
# when it seems that one has more "mother" than "father"; we can disregard this 
# filter for tumour. 
# 
# Note: recently Mutect2 has been extended to run multi-sample calling. 
# Now our approach would not be required anymore (all within Mutect2). 

# Take one sample, e.g., Set6_48
one_sample = sample_mutations$Set6_48

# Filters adopted
one_sample$FILTER %>% table %>% sort(decreasing = TRUE) %>% as_tibble()

# Cancel out "alleleBias"
one_sample$FILTER = gsub("alleleBias", "", one_sample$FILTER)

one_sample$FILTER %>% table %>% sort(decreasing = TRUE) %>% as_tibble()

one_sample = one_sample %>% 
  mutate(
    FILTER = ifelse(FILTER == "", "PASS", FILTER)
  )

one_sample$FILTER %>% table %>% sort(decreasing = TRUE) %>% as_tibble()

# Visualise the data
ggplot(one_sample, aes(VAF)) +
  geom_histogram(aes(fill = FILTER), binwidth = 0.01) +
  xlim(0, 1) +
  theme(legend.position = 'bottom') +
  facet_wrap(~FILTER)

# Most are PASS, retain only those
retain = "PASS"

one_sample = one_sample %>% filter(FILTER %in% !!retain)

ggplot(one_sample, aes(VAF)) +
  geom_histogram(aes(fill = FILTER), binwidth = 0.01) +
  xlim(0, 1) +
  theme(legend.position = 'bottom') +
  facet_wrap(~FILTER)

# Inspect the coverage and num_variants distribution (DP = depth, NV = number of reads with variant)
dp_plot = ggplot(one_sample, aes(DP)) +
  geom_histogram(aes(fill = FILTER), bins = 100) +
  theme(legend.position = 'bottom')  +
  scale_x_log10()

nv_plot = ggplot(one_sample, aes(NV)) +
  geom_histogram(aes(fill = FILTER), bins = 100) +
  theme(legend.position = 'bottom')  +
  scale_x_log10()

# scatter NV vs DP
dpnv_plot = ggplot(one_sample, aes(x = NV, y = DP)) +
  geom_point(aes(color = FILTER), size = .3) +
  geom_density2d() +
  theme(legend.position = 'bottom') +
  facet_wrap(~FILTER) +
  scale_y_log10()

# One figure
library(patchwork) # https://patchwork.data-imaginist.com/

(dp_plot/nv_plot) | dpnv_plot

# Take some filters - DP > 30, NV > 5, replot
one_sample = one_sample %>% 
  mutate(
    chop = DP < 30 | NV < 5
  )

dp_plot = ggplot(one_sample, aes(DP)) +
  geom_histogram(aes(fill = chop), bins = 100) +
  theme(legend.position = 'bottom')  +
  scale_x_log10()

nv_plot = ggplot(one_sample, aes(NV)) +
  geom_histogram(aes(fill = chop), bins = 100) +
  theme(legend.position = 'bottom')  +
  scale_x_log10()

# scatter NV vs DP
dpnv_plot = ggplot(one_sample, aes(x = NV, y = DP)) +
  geom_point(aes(color = chop), size = .3) +
  geom_density2d(
    data = one_sample %>% filter(!chop)
  ) +
  theme(legend.position = 'bottom') +
  scale_y_log10()

(dp_plot/nv_plot) | dpnv_plot

# How the data distribution changes for VAF
ggplot(
  one_sample,
  aes(VAF)
) +
  geom_histogram(aes(fill = chop), binwidth = 0.01) +
  xlim(0, 1) +
  theme(legend.position = 'bottom') +
  facet_wrap( ~ chop)

# VAF per chromosome - why the signal is different across chromosomes?
ggplot(
  one_sample %>% filter(!chop, VAF > 0.05),
  aes(VAF)
) +
  geom_histogram(binwidth = 0.01) +
  xlim(0, 1) +
  theme(legend.position = 'bottom') +
  facet_wrap( ~ chr)



