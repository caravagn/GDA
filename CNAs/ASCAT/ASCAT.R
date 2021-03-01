# Does not work it is not a package
# devtools::install_github('Crick-CancerGenomics/ascat')

devtools::build('ascat-master/ASCAT')
install.packages('./ascat-master/ASCAT_2.5.2.tar.gz', type = 'sources', ref = NULL)

library(ASCAT)

ascat.bc = ascat.loadData("ExampleData/Tumor_LogR.txt",
                          "ExampleData/Tumor_BAF.txt",
                          "ExampleData/Germline_LogR.txt", 
                          "ExampleData/Germline_BAF.txt")

plot_sample = function(x) {
  f1 = ascat.bc$Tumor_BAF[[x]] %>%
    as_tibble() %>%
    mutate(x = row_number()) %>%
    ggplot(aes(x = x, y = value, color = value)) +
    geom_point(size = .3) +
    scale_color_distiller(palette = 'Spectral') +
    theme_light()
  
  f2 = ascat.bc$Tumor_LogR[[x]] %>%
    as_tibble() %>%
    mutate(x = row_number()) %>%
    ggplot(aes(x = x, y = value, color = value)) +
    geom_point(size = .3) +
    scale_color_gradient2() +
    theme_light()
  
  ggpubr::ggarrange(f1, f2, ncol = 1, nrow = 2)
}

plot_sample("S1")

# Plot row data with ASCAT function
ascat.plotRawData(ascat.bc)

ascat.bc = ascat.aspcf(ascat.bc)

ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc) 

require(dplyr)

inp = read.table("Tumor_BAF.txt") %>% 
  filter(chrs == 1)

require(ggplot2)

inp_g = read.table("Germline_BAF.txt") %>% 
  filter(chrs == 1)

ggpubr::ggarrange(
  ggplot(inp, aes(x = pos, y = S1)) + geom_point(color = 'indianred3') + labs(title = 'Tumour'),
  ggplot(inp_g, aes(x = pos, y = S1)) + geom_point(color = 'steelblue') + labs(title = 'Germline'),
  nrow = 2
)

ggplot(inp %>% filter(S1 >0.1 & S1 <.9), aes(x = pos, y = ifelse(S1 < 0.5, .5+S1, S1))) + geom_point(color = 'indianred3') + labs(title = 'Tumour')



