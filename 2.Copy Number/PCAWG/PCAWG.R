require(dplyr)
require(ggplot2)

####### Load all the input data, per sample 

calls =
  lapply(
    list.files('samples/', full.names = TRUE),
    FUN = function(x) {
      s = gsub("samples//",
               '',
               gsub(".consensus.20170119.somatic.cna.annotated.txt", '', x))
      
      readr::read_tsv(x, col_types = readr::cols()) %>%
        mutate(sample_id = s)
    }
    ) %>%
  Reduce(f = bind_rows)

calls = calls %>% 
  select(chromosome, start, end, total_cn, major_cn, minor_cn, sample_id) %>% 
  distinct(chromosome, start, end, .keep_all= TRUE)

####### Load this table with the association between mutations and tumour type
drv = readr::read_tsv("TableS3_panorama_driver_mutations_ICGC_samples.public.tsv") %>% 
  select(sample_id, ttype)

calls = calls %>% left_join(drv) %>% filter(!is.na(ttype)) 


####### Filter out only clonal CNA segments using as reference the total_cn column,
# remove cray high CN
calls = calls %>% 
  filter(total_cn < 10)

# These are the tumour types we are considering 
calls$ttype %>% unique()

# This is what I created drafting this simple code
# calls %>% readr::write_csv("CNA_pcawg.csv")

example_chr = '1'

# Play around with this dataset and try to make the following visualisation
# - total number of segments per tumour type
# - per chromosome barplot of total CN per tumour type
# - barplot as above, but retaining only diploid segments (1 = Major, 1 = minor)
# - overall (whole-genome) plot of 
# - rank chromosomes by number or size of LOH events (0 = minor) 

calls %>% 
  filter(chromosome == example_chr) %>% 
  ggplot(aes(total_cn, fill = ttype)) +
  geom_bar() +
  scale_fill_manual(values = pals::cols25())

calls %>% 
  filter(!(major_cn == 1 && minor_cn == 1), chromosome == example_chr) %>% 
  ggplot(aes(interaction(major_cn, minor_cn), fill = ttype)) +
  geom_bar() +
  scale_fill_manual(values = pals::cols25()) +
  facet_wrap(~chromosome) 

calls %>% 
  filter(!(major_cn == 1 && minor_cn == 1)) %>% 
  ggplot(aes(total_cn, fill = ttype)) +
  geom_bar() +
  scale_fill_manual(values = pals::cols25()) +
  facet_wrap(~chromosome) +
  scale_x_discrete()

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


  
  
  