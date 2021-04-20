# load required libraries and sources
library("SparseSignatures")

# load the data
load("background.RData")
load("trinucleotides_counts.RData")

# settings
K = 2:7
nmf_runs = 10
my_seed_starting_beta = 789

# fit the initial betas for each configuration
initial_betas = startingBetaEstimation(x=trinucleotides_counts,K=K,background_signature=background,normalize_counts=FALSE,nmf_runs=nmf_runs,seed=my_seed_starting_beta,verbose=TRUE)
save(initial_betas,file="results/initial_betas.RData")
