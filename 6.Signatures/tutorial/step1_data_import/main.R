# load the required libraries
library("data.table")
library("BSgenome.Hsapiens.1000genomes.hs37d5")
source("data.import.R")

# load data
load(file="raw_data.RData")

# make trinucleotide counts matrix
trinucleotides_counts = get.SBS.counts(data=raw_data,reference=BSgenome.Hsapiens.1000genomes.hs37d5)

# save the results
save(trinucleotides_counts,file="trinucleotides_counts.RData")
