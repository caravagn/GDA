# load required libraries
library("SparseSignatures")
library("data.table")
library("ggplot2")
library("gridExtra")
library("lsa")
"signatures.plot.v2" = function (beta, useColNames = TRUE, mutation_categories = NULL, firstBackground = TRUE, xlabels = TRUE) {
    if (firstBackground) {
        rownames(beta) <- c("Background", paste0("Signature ", 
            1:(nrow(beta) - 1)))
    }
    else {
        ###rownames(beta) <- paste0("Signature ", 1:nrow(beta))
    }
    if (!useColNames) {
        colnames(beta) <- sort(mutation_categories$cat)
    }
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    x[, `:=`(Context, paste0(substr(cat,1,1), ".", substr(cat, 
        7, 7)))]
    x[, `:=`(alt, paste0(substr(cat,3,3), ">", substr(cat, 
        5, 5)))]
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = Context, 
            y = value, fill = alt), stat = "identity", position = "identity") + 
            facet_wrap(~alt, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") + 
            ylab("Frequency of mutations")
        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        }
        glist[[i]] <- plt
    }
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
}

# read the results
load("trinucleotides_counts.RData")
load("results_inference/nmf_brunet.RData")
load("results_inference/SigProfiler.RData")
load("results_inference/SparseSignatures.RData")

# read the signatures by COSMICv3.1
load("COSMIC_v3.1_June2020.RData")
cosmic_signatures = COSMIC_v3.1_June2020[,colnames(trinucleotides_counts)]

# make plots of optimal number of signatures
RANGE = 2:7
VALUE = as.numeric(nmf_brunet$goodness_fit[,"silhouette.consensus"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("NMF (Silhouette Consensus)") + xlab("Number of Signatures") + ylab("Mean Silhouette Consensus (10 restarts)")

RANGE = 2:7
VALUE = as.numeric(nmf_brunet$goodness_fit[,"evar"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("NMF (Explained Variance)") + xlab("Number of Signatures") + ylab("Mean Explained Variance (10 restarts)")

RANGE = 2:7
VALUE = as.numeric(SigProfiler$rank[,"Stability"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SigProfiler (Stability)") + xlab("Number of Signatures") + ylab("Mean Stability (10 restarts)")

RANGE = 2:7
VALUE = as.numeric(SigProfiler$rank[,"MeanL2"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SigProfiler (L2 Norm)") + xlab("Number of Signatures") + ylab("Mean L2 Norm (10 restarts)")

RANGE = 2:7
VALUE = as.numeric(SparseSignatures$rank[,"Median cross-validation error"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SparseSignatures (Cross Validation Error)") + xlab("Number of Signatures") + ylab("Median Cross Validation Error (10 restarts)")

RANGE = 2:7
VALUE = as.numeric(SparseSignatures$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SparseSignatures (Goodness of fit)") + xlab("Number of Signatures") + ylab("Goodness of fit (percentage of predictions with >0.95 cosine similarity)")

# NMF
alpha = nmf_brunet$alpha["Rank 5",][[1]]
beta = nmf_brunet$beta["Rank 5",][[1]]
silhouette_consensus = as.numeric(nmf_brunet$goodness_fit[,"silhouette.consensus"])
names(silhouette_consensus) = paste0("Signature ",2:7)
explained_variance = as.numeric(nmf_brunet$goodness_fit[,"evar"])
names(explained_variance) = paste0("Signature ",2:7)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,2,5,1,4)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
names = c("NMF1 (SBS2 - 0.91)","NMF2 (SBS3 - 0.94)","NMF3 (SBS6 - 0.70)","NMF4 (SBS13 - 0.77)","NMF5 (SBS44 - 0.84)")
colnames(alpha) = names
rownames(beta) = names
NMF = list()
NMF[["alpha"]] = alpha
NMF[["beta"]] = beta
NMF[["Silhouette_Consensus"]] = silhouette_consensus
NMF[["Explained_Variance"]] = explained_variance

# SigProfiler
alpha = SigProfiler$alpha["Rank 5",][[1]]
beta = SigProfiler$beta["Rank 5",][[1]]
Stability = as.numeric(SigProfiler$rank[,"Stability"])
names(Stability) = paste0("Signature ",2:7)
MeanL2 = as.numeric(SigProfiler$rank[,"MeanL2"])
names(MeanL2) = paste0("Signature ",2:7)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(5,2,1,4,3)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
names = c("SIP1 (SBS2 - 0.90)","SIP2 (SBS3 - 0.96)","SIP3 (SBS13 - 0.78)","SIP4 (SBS18 - 0.68)","SIP5 (SBS44 - 0.82)")
colnames(alpha) = names
rownames(beta) = names
SigProfiler = list()
SigProfiler[["alpha"]] = alpha
SigProfiler[["beta"]] = beta
SigProfiler[["Stability"]] = Stability
SigProfiler[["Mean_L2"]] = MeanL2

# SparseSignatures
alpha = SparseSignatures$alpha["Rank 5",][[1]]
beta = SparseSignatures$beta["Rank 5",][[1]]
median_cv_error = as.numeric(SparseSignatures$rank[,"Median cross-validation error"])
names(median_cv_error) = paste0("Signature ",2:7)
goodness_fit = as.numeric(SparseSignatures$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
names(goodness_fit) = paste0("Signature ",2:7)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,2,1,5,4)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
names = c("SPS1 (SBS2 - 0.77)","SPS2 (SBS3 - 0.88)","SPS3 (SBS5 - 0.98)","SPS4 (SBS13 - 0.78)","SPS5 (SBS4 - 0.81)")
colnames(alpha) = names
rownames(beta) = names
SparseSignatures = list()
SparseSignatures[["alpha"]] = alpha
SparseSignatures[["beta"]] = beta
SparseSignatures[["Median_CV_Error"]] = median_cv_error
SparseSignatures[["Goodness_Fit"]] = goodness_fit

# save results
save(NMF,file="final_results/NMF.RData")
save(SigProfiler,file="final_results/SigProfiler.RData")
save(SparseSignatures,file="final_results/SparseSignatures.RData")

# make plots of discovered signatures
signatures.plot.v2(NMF$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(SigProfiler$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(SparseSignatures$beta,firstBackground=FALSE,xlabels=FALSE)
