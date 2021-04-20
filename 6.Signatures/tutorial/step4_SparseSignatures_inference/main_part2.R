# load required libraries and sources
library("SparseSignatures")

# load the data
load("background.RData")
load("trinucleotides_counts.RData")
load("results/initial_betas.RData")

# settings
K = 2:7
lambda_values_alpha = c(0.000)
lambda_values_beta = c(0.000,0.025)
cross_validation_entries = 0.01
cross_validation_repetitions = 10
num_processes = 10
my_seed_cv = 987
log_file = ""

# performing cross-validation for a grid search to estimate the optimal number of signatures
cross_validation = nmfLassoCV(x=trinucleotides_counts,K=K,starting_beta=initial_betas,normalize_counts=FALSE,lambda_values_alpha=lambda_values_alpha,lambda_values_beta=lambda_values_beta,cross_validation_entries=cross_validation_entries,cross_validation_repetitions=cross_validation_repetitions,num_processes=num_processes,seed=my_seed_cv,verbose=TRUE,log_file=log_file)
save(cross_validation,file="results/cross_validation.RData")
