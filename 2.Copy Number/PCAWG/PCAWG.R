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

drv = readr::read_tsv("~/Downloads/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv") %>% 
  select(sample_id, ttype)

calls = calls %>% left_join(drv) %>% filter(!is.na(ttype)) 

calls = calls %>% filter(total_cn < 10)

calls$ttype %>% unique()

# ttypes = calls %>% distinct(sample_id, ttype) %>% group_by(ttype) %>% summarise(n = n()) %>% arrange(desc(n)) %>%
#   pull(ttype)

calls %>% readr::write_csv("CNA_pcawg.csv")

example_chr = '1'

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
