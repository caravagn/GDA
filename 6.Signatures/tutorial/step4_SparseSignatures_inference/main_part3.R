# data structure to save final results
cross_validation_summary = array(NA,c(6,2))
rownames(cross_validation_summary) = paste0("Rank ",2:7)
colnames(cross_validation_summary) = c("Median_CV_Error","Best_Lambda")

# save results
load(file="results/cross_validation.RData")
for(i in 1:6) {
    cv_medians = NULL
    for(j in names(cross_validation$grid_search_mse[,,i])) {
        cv_medians = c(cv_medians,median(cross_validation$grid_search_mse[,j,i][[1]]))
    }
    cross_validation_summary[i,"Median_CV_Error"] = min(cv_medians)
    cross_validation_summary[i,"Best_Lambda"] = as.numeric(gsub("_Lambda_Beta","",names(cross_validation$grid_search_mse[,,i])[which(cv_medians==min(cv_medians))[1]]))
}
save(cross_validation_summary,file="results/cross_validation_summary.RData")
