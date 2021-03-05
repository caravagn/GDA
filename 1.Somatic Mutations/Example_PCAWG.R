# Download the PCAWG cohort data
#
# -> https://dcc.icgc.org/releases/PCAWG/
# -> https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
# -> https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz
#
# ~2700 WGS pan-cancer, ~1GB of somatic mutations calls.
#
# # Samples for which we will have CNA data (next lecture)
# example_samples = readr::read_csv("2.Copy Number/PCAWG/CNA_pcawg.csv") %>%
#   pull(sample_id) %>%
#   unique
# #
# # # I prepare data for this lecture (after download)
# muts = readr::read_tsv("~/Downloads/final_consensus_passonly.snv_mnv_indel.icgc.public.maf")
# muts_example_samples = muts %>%
#   filter(Tumor_Sample_Barcode %in% example_samples) %>%
#   select(
#     Chromosome,
#     Start_position,
#     End_position,
#     Reference_Allele,
#     Tumor_Seq_Allele1,
#     t_alt_count,
#     t_ref_count,
#     Tumor_Sample_Barcode,
#     everything()
#   )
#
# # Give better names for next lecture
# colnames(muts_example_samples)[1:8] = c('chr', 'from', 'to', 'ref', 'alt', 'NV', 'NR', 'sample')
#
# muts_example_samples = muts_example_samples %>%
#   mutate(
#     DP = NV + NR,
#     VAF = NV / DP,
#     chr = paste0('chr', chr),
#     alt = ifelse(ref == alt, Tumor_Seq_Allele2, alt)
#   ) %>%
#   select(chr, from, to, ref, alt, DP, NV, VAF, sample, everything())
#
#
# lapply(example_samples,
#        function(x) {
#          muts_example_samples %>%
#            filter(sample == x) %>%
#            readr::write_csv(file = paste0('1.Somatic Mutations/PCAWG/samples/', x, '.csv'))
#        })
#
# ttype %>% readr::write_csv('1.Somatic Mutations/PCAWG/Types.csv')

#
# Lecture
#

# Themes for ggplot2
install.packages("hrbrthemes")

ggplot2::reset_theme_settings()
ggplot2::theme_set(hrbrthemes::theme_ipsum(base_size = 8) +
                     theme(legend.position = 'bottom'))

# PCAWG, 27 samples
all_files = list.files('1.Somatic Mutations/PCAWG/samples/', full.names = TRUE)

all_mutations = lapply(all_files, readr::read_csv)
all_mutations = lapply(all_mutations,
                       function(x) {
                         x$i_NumCallers = paste0(x$i_NumCallers)
                         x$i_VAF = paste0(x$i_VAF)
                         x
                       })

all_mutations = Reduce(bind_rows, all_mutations)

# ~500.000 mutations, add tumor type
ttype = readr::read_csv('1.Somatic Mutations/PCAWG/Types.csv')

all_mutations = all_mutations %>% full_join(ttype)

# Set mutation type
all_mutations = all_mutations %>%
  mutate(mtype = ifelse(
    ref %in% c("A", "C", "T", "G") & alt %in% c("A", "C", "T", "G"),
    "SNV",
    "indel"
  ))

# Distribution of SNVs vs indels across ttype
ggplot(
  all_mutations %>%
    group_by(mtype, ttype) %>%
    summarise(n = n()),
  aes(x = '', y = n, fill = mtype)
) +
  geom_bar(stat = 'identity', position = position_fill()) +
  facet_wrap( ~ ttype) +
  ggsci::scale_fill_lancet() +
  coord_polar("y", start = 0) +
  theme(axis.text.x = element_blank())

# x = all_mutations
# all_mutations = all_mutations %>% filter(ttype == "Breast-AdenoCa" )

# Mutational burden per sample
order_ttype = all_mutations  %>% 
  # filter(ttype == "Breast-AdenoCa" ) %>% 
  group_by(sample) %>%
  summarise(n = median(n(), na.rm = TRUE), ttype = ttype[1]) %>% 
  group_by(ttype) %>%
  summarise(n = median(n, na.rm = TRUE)) %>% 
  arrange(n)

all_mutations  %>% 
  # filter(ttype == "Breast-AdenoCa" ) %>% 
  group_by(sample, ttype) %>%
  summarise(n = n(), DP = median(DP, na.rm = TRUE)) %>% 
  ggplot(aes(x = ttype, y = n, size = DP)) +
  geom_jitter() +
  scale_y_log10() +
  scale_size_continuous() +
  scale_x_discrete(limits = order_ttype$ttype) +
  geom_point(
    data = order_ttype,
    aes(x = ttype, y = n),
    inherit.aes = FALSE,
    color = 'indianred3',
    size = 3,
    pch = 17
  ) +
  theme(
    axis.text.x = element_text(angle = 90)
  )

# Reproduce this plot: 
# 
# https://github.com/cBioPortal/icebox/issues/78 
# 
# using the full PCAWG cohort
