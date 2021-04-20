# load required libraries and sources
library("nnls")

# load data and results
load("trinucleotides_counts.RData")

# set the seed
set.seed(654)

# fit final results
beta = array(list(),c(6,1))
rownames(beta) = paste0("Rank ",2:7)
colnames(beta) = "Results"
alpha = beta
for(k in 2:7) {
    res = read.table(paste0("results/SBS96/All_solutions/SBS96_",k,"_Signatures/SBS96_S",k,"_Signatures.txt"),header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
    curr_alpha = array(NA,c(nrow(trinucleotides_counts),k))
    rownames(curr_alpha) = rownames(trinucleotides_counts)
    colnames(curr_alpha) = paste0("Signature ",1:ncol(curr_alpha))
    curr_beta = array(NA,c(k,ncol(trinucleotides_counts)))
    rownames(curr_beta) = colnames(curr_alpha)
    colnames(curr_beta) = colnames(trinucleotides_counts)
    for(i in 2:ncol(res)) {
        for(j in 1:nrow(res)) {
            curr_beta[(i-1),res[j,"MutationType"]] = res[j,i]
        }
    }
    for(j in 1:nrow(curr_alpha)) {
        curr_alpha[j,] = nnls(t(curr_beta),as.vector(trinucleotides_counts[j,]))$x
    }
    alpha[[(k-1),1]] = curr_alpha
    beta[[(k-1),1]] = curr_beta
}
rank = read.csv("results/SBS96/All_solutions_stat.csv")
rank = rank[,c("Stability..Avg.Silhouette.","Mean.Sample.L2.")]
Stability = as.numeric(rank$Stability..Avg.Silhouette.)
MeanL2 = as.numeric(gsub("%","",rank$Mean.Sample.L2.))/100
rank = cbind(Stability,MeanL2)
rownames(rank) = 2:7

# save results
SigProfiler = list(alpha=alpha,beta=beta,rank=rank)
save(SigProfiler,file="SigProfiler.RData")
