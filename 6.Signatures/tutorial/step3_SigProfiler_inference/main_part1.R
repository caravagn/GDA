# load required libraries and sources
library("SigProfilerExtractorR")

# save input data
load("trinucleotides_counts.RData")
data = cbind(t(t(colnames(trinucleotides_counts))),t(trinucleotides_counts))
colnames(data)[1] = "Mutation Types"
write.table(data,file="input_data/data.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

# set the seed
set.seed(456)

# perform inference
sigprofilerextractor(input_type="table",output="results",input_data="input_data/data.txt",minimum_signatures=2,maximum_signatures=7,nmf_replicates=10,init="random",exome=FALSE,cpu=10)
