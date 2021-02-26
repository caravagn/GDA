devtools::install_github('Crick-CancerGenomics/ascat')
devtools::build()
install.packages('./ascat-master/ASCAT_2.5.2.tar.gz', type = 'sources', ref = NULL)


library(ASCAT)

ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt","Germline_LogR.txt", "Germline_BAF.txt")

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



