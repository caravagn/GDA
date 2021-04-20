# load required libraries and sources
library("SparseSignatures")
library("nnls")
library("lsa")

# load the data
load("background.RData")
load("trinucleotides_counts.RData")
load("results/initial_betas.RData")
load("results/cross_validation_summary.RData")

# set the seed
set.seed(789)

# get final results
K = 2:7
beta = array(list(),c(length(K),1))
rownames(beta) = paste0("Rank ",K)
colnames(beta) = "Results"
alpha = beta
goodness_fit_mean = NULL
goodness_fit_percent = NULL
for(k in 2:7) {
    # perform fit for beta
    res = nmfLasso(x=trinucleotides_counts,K=k,beta=initial_betas[[(k-1),1]],normalize_counts=FALSE,lambda_rate_alpha=0.000,lambda_rate_beta=cross_validation_summary[(k-1),"Best_Lambda"],seed=k,verbose=FALSE)
    curr_beta = res$beta
    rownames(curr_beta) = paste0("Signature ",1:nrow(curr_beta))
    colnames(curr_beta) = colnames(trinucleotides_counts)
    # perform fit for alpha
    curr_alpha = array(NA,c(nrow(trinucleotides_counts),k))
    rownames(curr_alpha) = rownames(trinucleotides_counts)
    colnames(curr_alpha) = paste0("Signature ",1:ncol(curr_alpha))
    for(i in 1:nrow(trinucleotides_counts)) {
        curr_alpha[i,] <- nnls(t(curr_beta),as.vector(trinucleotides_counts[i,]))$x
    }
    alpha[[(k-1),"Results"]] = curr_alpha
    beta[[(k-1),"Results"]] = curr_beta
    # estimate goodness of fit
    predicted_counts = curr_alpha%*%curr_beta
    curr_goodness_fit = NULL
    for(i in 1:nrow(trinucleotides_counts)) {
        curr_goodness_fit = c(curr_goodness_fit,as.numeric(cosine(as.numeric(trinucleotides_counts[i,]),as.numeric(predicted_counts[i,]))))
    }
    goodness_fit_mean = c(goodness_fit_mean,mean(curr_goodness_fit))
    goodness_fit_percent = c(goodness_fit_percent,(length(which(curr_goodness_fit>0.95))/length(curr_goodness_fit)))
    cat((k-1)/6,"\n")
}
rank_estimation = cbind(cross_validation_summary[,"Median_CV_Error"],goodness_fit_mean,goodness_fit_percent)
colnames(rank_estimation) = c("Median cross-validation error","Goodness of fit (mean cosine similarity)","Goodness of fit (percent >0.95 cosine similarity)")

# save results
SparseSignatures = list(alpha=alpha,beta=beta,rank=rank_estimation)
save(SparseSignatures,file="SparseSignatures.RData")
