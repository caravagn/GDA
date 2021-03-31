# AT HOME
# 
# Download the PCAWG cohort data
#
# -> https://dcc.icgc.org/releases/PCAWG/
# -> https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
# -> https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz
#
# ~2700 WGS pan-cancer, ~1GB of somatic mutations calls.

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

##########################################################################################
# Lecture - Somatic Mutation Calling
##########################################################################################

# I assembled some PCAWG cases for you (somatic calls)
PCAWG_example = "https://www.dropbox.com/s/stkkocdjymmtcp8/PCAWG_somatic_mutations.zip?raw=1"

# Download, load and cancel data
download.file(PCAWG_example, "PCAWG_somatic_mutations.zip")

unzip("PCAWG_somatic_mutations.zip")

file.remove("PCAWG_somatic_mutations.zip")

########################################################################################## 
# Download, unzip and cancel data, then load all the data in a tibble. 
# There should be ~500.000 mutations, add tumor type joining data file Types.csv, 
# and add a categorical flag "SNV", "indel" based on the type of mutations identified 
# (use reference and alternative alleles to determine that)
##########################################################################################
all_files = list.files('./samples/', full.names = TRUE)

all_mutations = lapply(all_files, readr::read_csv)
all_mutations = lapply(all_mutations,
                       function(x) {
                         x$i_NumCallers = paste0(x$i_NumCallers)
                         x$i_VAF = paste0(x$i_VAF)
                         x
                       })

all_mutations = Reduce(bind_rows, all_mutations)

##########################################################################################
# Visualise the distribution of SNVs/indels across tumour types, with a piechart/ barplot
##########################################################################################
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

##########################################################################################
# Compute mutational burden per sample, split by type of mutation
##########################################################################################
all_mutations %>% 
  group_by(sample, mtype) %>%
  summarise(n = n(), ttype = ttype[1]) %>% 
  arrange(desc(n)) %>% 
  ggplot()+
  geom_bar(aes(x = sample, fill = mtype, y = n), stat = 'identity') +
  # scale_y_log10() +
  coord_flip()

##########################################################################################
# Define the substitutions for an SNV as C>T, compute their distribution per cancer type.
# Which tumour has more mutations? Can you guess why?
##########################################################################################
all_mutations  %>% 
  filter(mtype == 'SNV') %>% 
  mutate(substitution = paste0(ref, '>', alt)) %>% 
  group_by(ttype, substitution) %>%
  summarise(n = median(n(), na.rm = TRUE), ttype = ttype[1]) %>% 
  arrange(desc(n)) %>% 
  ggplot() +
  geom_bar(aes(x = ttype, n, fill = substitution), stat = 'identity') +
  coord_flip()
  
##########################################################################################
# Compute the mutational burden per tumour type, sort types by that, and use a scatter
# where you can annotate each sample by median depth of sequencing
##########################################################################################

order_ttype = all_mutations  %>% 
  # filter(ttype == "Breast-AdenoCa" ) %>% 
  group_by(sample) %>%
  summarise(n = median(n(), na.rm = TRUE), ttype = ttype[1]) %>% 
  group_by(ttype) %>%
  summarise(n = median(n, na.rm = TRUE)) %>% 
  arrange(n)

all_mutations  %>% 
  group_by(sample, ttype) %>%
  summarise(n = n(), DP = median(DP, na.rm = TRUE)) %>% 
  ggplot(aes(x = ttype, y = n, size = DP)) +
  geom_jitter() +
  # scale_y_log10() +
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


##########################################################################################
# Compute which genes has the highest number of "Missense_Mutation" flags. Do you
# see any gene you ever heard of?
##########################################################################################
all_mutations %>% 
  filter(Variant_Classification == "Missense_Mutation") %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))


##########################################################################################
# Reproduce this plot: 
# 
# https://github.com/cBioPortal/icebox/issues/78 
##########################################################################################

##########################################################################################
# Remake these plots with all the PCAWG cohort. Note that if you download the data you
# will have to process data columns etc (I polished those for you now).
##########################################################################################

