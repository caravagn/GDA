library(tidyverse)

setwd("10.REVOLVER/")

############ ############ ############ ############ ############ ############  
# Part 1 - playing with REVOLVER
############ ############ ############ ############ ############ ############ 
# devtools::install_github('caravagnalab/revolver')
library(revolver)

# Get data from the paper: install.github("caravagn/evoverse.datasets")
data('TRACERx_NEJM_2017_REVOLVER', package = 'evoverse.datasets')

TRACERx_NEJM_2017_REVOLVER

TRACERx_NEJM_2017_REVOLVER %>% plot

TRACERx_NEJM_2017_REVOLVER$dataset$CRUK0002 %>% View

# Refit the cohort
options(easypar.parallel = FALSE)

refit = revolver_fit(TRACERx_NEJM_2017_REVOLVER)

refit = revolver_cluster(refit)

# Inspect some data plot
revolver:::plot_patient_trees(refit, "CRUK0002")
revolver:::plot_patient_data(refit, "CRUK0002")
revolver:::plot_patient_trees_scores(refit, "CRUK0002")

# Focus on WEX vs WGS
revolver::plot_patient_CCF_histogram(refit, patient = 'CRUK0002')

# What is clonal vs subclonal
revolver::plot_drivers_clonality(refit)

# Trajectories per cluster
revolver:::plot_trajectories_per_cluster(refit, min_counts = 2)

# Plot cluster
revolver::plot_clusters(refit, cutoff_drivers = 3)

# Distance as dendrogram
revolver::plot_dendrogram(refit)

############ ############ ############ ############ ############ ############ 
# Part 2 - cluster enrichment and survival analysis 
############ ############ ############ ############ ############ ############ 

# # Revolver input dataset
# breast_fit = revolver_cohort(evoverse.datasets::YATES_BREAST_NATMED_2015 %>% 
#                                filter(!(patientID %in% c('PD14779a'))), 
#                              MIN.CLUSTER.SIZE = 0, 
#                              annotation = "Yates breast cancer data (multi-region)")
# 
# # breast_fit = remove_patients(breast_fit, 'PD14779a')
# 
# # Build trees and fit the data
# breast_fit = breast_fit %>% compute_clone_trees(sspace.cutoff = 1000, 
#                                                 n.sampling = 500, 
#                                                 store.max = 100)
# breast_fit = breast_fit %>% revolver_fit()
# 
# # Cluster
# breast_fit = breast_fit %>% revolver_cluster()

load('Breast-cancer-survival.RData', verbose = TRUE)
names(test.set)

# Inspect the data
tibble::as_tibble(test.set$data)
tibble::as_tibble(test.set$clinical_data)

# Classification of the test set
classify = function(patient)
{
  group = 'C1'
  
  if(patient['TP53'] == 1) 
  {
    if(patient['+8q'] == 1)
    { 
      if ((patient['+1q'] == 1) || (patient['-17p'] == 1) ) group = 'C5'
      else group = 'C3'
    } else group = 'C2'
  } else { 
    
    if(patient['-16q'] == 1) 
    {
      if(patient['+8q'] == 1)
      {
        if(patient['PIK3CA'] == 1) group = 'C6'
        else group = 'C4'
      }
    }
  }  
  
  return(group)
}


# Then we begin by classifying our samples
data = test.set$data
clinical_data = test.set$clinical_data
classification = rep(NA, nrow(data))
names(classification) = rownames(data)
for(i in 1:nrow(data)) classification[i] = classify(data[i, ])
print(head(classification, 15))

# We use ``pheatmap`` to plot the data, so we create annotations for 
# IntClust and Luminal-Basal classification, and define some colours for the plots 

# The pheatmap annotation data frame
annotations = as.data.frame(classification, stringsAsFactors = FALSE)
colnames(annotations) = 'cluster'
head(annotations)

# Colours of the plot
colors = c("#E41A1C","#FFFF33", "#984EA3", "#377EB8", "#FF7F00" , "#4DAF4A")
names(colors) = c("C1", "C2", "C3", "C4", "C5","C6")

# Piechart of mixing proportions
slices = as.vector(table(classification))
lbls = as.vector(names(table(classification)))
pct = round(slices/sum(slices)*100)
lbls =  paste(lbls, pct)
lbls =  paste(lbls, "%", sep="")
pie(slices, labels = lbls, col = colors)
legend("topright", lbls, cex=0.8, fill=colors)

# We annotate the status of the other classifications as well (IntClust and PAM50) on our test set.

annotations = cbind(annotations, 
                    IntClust = c(clinical_data[rownames(annotations), 'INTCLUST']))

# With not much fantasy we use the term Classification to denote the PAM50 grouping of breast cancers
annotations = cbind(annotations, 
                    Classification = c(clinical_data[rownames(annotations), 'PAM50']))

# We control for some NAs and factors
annotations$IntClust = addNA(annotations$IntClust)
annotations$IntClust = factor(
  annotations$IntClust,
  levels = c("1", "2" ,"3","4","4ER-","4ER+","5","6","7","8","9","10" ,"unknown")
)

# And annotate samples with unknown classification
id.na = which(is.na(annotations$IntClust))
annotations[id.na,'IntClust'] = "unknown"

# We define colors for annotations

# IntClust
bpalette = RColorBrewer::brewer.pal(8, 'Accent')
ncols = length(unique(annotations$IntClust))
col.IntClust = colorRampPalette(bpalette)(ncols)
names(col.IntClust) = c("1", "2" ,"3","4","4ER-","4ER+","5","6","7","8","9","10" ,"unknown")    

# PAM50
col.Classification = RColorBrewer::brewer.pal(3, 'Dark2')[1:3]
names(col.Classification) = c("Luminal A", "Luminal B", "Basal")

# We create a list of colors for ``pheatmap``.
annotation_colors = list(
  cluster = colors, 
  Classification = col.Classification, 
  IntClust = col.IntClust)

# and we order the columns according to the clusters
annotations = annotations[order(annotations$cluster),  , drop = FALSE]
annotations = annotations[
  with(annotations, order(annotations$cluster, annotations$Classification)), ]

annotations = annotations[
  with(annotations, order(annotations$Classification, annotations$IntClust)), ]

annotations.ordered = NULL
for(c in names(colors))
{
  sd = annotations[annotations$cluster == c, , drop = FALSE ]
  annotations.ordered = rbind(annotations.ordered, sd)
}
print(head(annotations))

# Because we have a large number of samples, we transpose the data.
tdata = t(data)
tdata = tdata[, rownames(annotations.ordered)]

pheatmap::pheatmap(
  tdata, 
  cellheight = 13, 
  cellwidth = 0.3, 
  legend = FALSE,
  main = "Classification of METABRIC/BRCA2012 LumA-LumB-Basal samples", 
  cluster_rows = F,
  cluster_cols = F,
  annotation_col = annotations.ordered,
  color = c('gray88', 'midnightblue'),
  annotation_colors = annotation_colors,
  show_colnames = F)

# Enrichment analysis - via Fisher's Exact Test.
classification = as.data.frame(classification, stringsAsFactors = FALSE)

# Enrichment IntClust
data.test = annotations[(annotations$cluster!="unknown"),c("IntClust", "cluster")]
head(data.test)

# Contingency table
ctable = function(G, Clust)
{
  table = table(data.test)
  diff = setdiff(colnames(table), G)
  diff2 = setdiff(rownames(table), Clust)
  ct = matrix(0, nrow = 2, ncol = 2)
  rownames(ct) = c(Clust, paste("Not", Clust))
  colnames(ct) = c(G, paste("Not", G))
  
  ct[1, 1] = table[Clust, G]
  ct[1, 2] =  sum(table[Clust, diff])
  ct[2, 1] =  sum(table[diff2,G])
  ct[2, 2] =  sum(table[diff2,diff])
  print(sum(ct))
  
  return(ct)  
}

IntClust = c("1", "2" , "3","4","4ER-","4ER+","5","6","7","8","9","10")
groups = c("C1", "C2", "C3",  "C4",  "C5", "C6")

# Matrix of p-values
p.value = matrix(0, length(groups), length(IntClust))
rownames(p.value) = groups
colnames(p.value) = IntClust

# Matrix of odds.ratios
odds.ratio = matrix(0, length(groups), length(IntClust))
rownames(odds.ratio) = groups
colnames(odds.ratio) = IntClust

for (i in groups){
  for(j in IntClust){
    cat('\nTesting ', i, '\n')
    print(ctable(i,j))
    
    # Fisher test
    test = fisher.test(ctable(i,j), alternative = "greater")
    
    p.value[i,j] = test$p.value
    odds.ratio[i, j] = test$estimate
  }
}

p = as.matrix(p.value)
o = round(odds.ratio, 3)


# The results of the test can be shown, and we can apply Bonferroni correction for $t=72$ tests.
print(p)

# Apply Bonferroni correction
p < (0.05/72)

# Correlation plot
corrplot::corrplot(p.value,
         p.mat = p.value,
         is.corr = FALSE,
         sig.level = -1,
         insig = "p-value", 
         cl.lim = c(0,1), 
         tl.cex = 1, 
         cl.cex = 1, 
         tl.col = "black")
`

# And then we do the same thing for the enrichment of Pam50 subtypes.

data.test = annotations[,c("Classification", "cluster")]
Class = c("Luminal A", "Luminal B", "Basal")

p.value = matrix(0, length(groups), length(Class))
rownames(p.value) = groups
colnames(p.value) = Class

odds.ratio = matrix(0, length(groups), length(Class))
rownames(odds.ratio) = groups
colnames(odds.ratio) = Class

for (i in groups) {
  for(j in Class) {
    cat('\nTesting ', i, '\n')
    print(ctable(i,j))
    
    test = fisher.test(ctable(i,j), alternative = "greater")
    
    p.value[i,j] = test$p.value
    odds.ratio[i, j] = test$estimate
  }
}

p = as.matrix(p.value)
o = round(odds.ratio, 3)

# And then correct for multiple testing as before, with $t=18$ tests. Bonferroni correction
p < (0.05/18)

corrplot::corrplot(p.value,
         p.mat = p.value,
         is.corr = FALSE,
         sig.level = -1,
         insig = "p-value", 
         cl.lim = c(0,1), 
         tl.cex = 2, 
         cl.cex = 1, 
         tl.col = "black",
         mar = c(0,0,1,0))


# Survival Analysis - We carry this out with two standard packages.

library(survival)

# install.packages('survminer')
library(survminer)

# We use overall survival as readout
clinical_data = clinical_data[,c('OS_MONTHS', 'OS_STATUS')]   
clinical_data = clinical_data[complete.cases(clinical_data),] 

clinical_data$OS_STATUS = sapply(
  clinical_data$OS_STATUS, function(x) if(x=="DECEASED") 1 else 0)

# We take as colours the one used in the main figures of the paper, and we start creating the structures that we need to use Cox regression.

groups =  c("C1", "C2", "C3", "C4", "C5", "C6")
colors = c("#E41A1C","#FFFF33", "#984EA3", "#377EB8", "#FF7F00" , "#4DAF4A")
names(colors) = groups

# We can analyze the samples for which we have all the required data
patients = intersect(rownames(clinical_data),  rownames(classification)) 

# We use a copy of clinic_data
cclinic_data = clinical_data[patients, , drop = FALSE]
classification = classification[patients, , drop = FALSE]
cclinic_data = cbind(cclinic_data, cluster = classification) 
colnames(cclinic_data) = c( "time", "death", "cluster")

# These will be used for survival computation
cclinic_data$time = as.numeric(cclinic_data$time)
cclinic_data$death = as.numeric(cclinic_data$death)
cclinic_data = cclinic_data[order(cclinic_data$cluster),]

# We begin the actual analysis of survival. First we plot the survival curve
# of each group separately.

for(c in groups)
{
  cd = cclinic_data[cclinic_data$cluster == c, ] 
  fit = survfit(Surv(cd$time, cd$death)~cd$cluster)
  
  pl = ggsurvplot(fit,
                  censor=F,
                  data = cd,  
                  surv.median.line = "hv",
                  linetype = 1, 
                  pval = TRUE,
                  ylab = 'Censored obs.',
                  ggtheme = theme_light(),
                  title = paste('Cluster', c),
                  color = colors[c],
                  ncensor.plot = TRUE)
  print(pl)
}

# Then we perform an overall survival analysis and compute a global p-value.

# Survival for all groups
fit = survfit(Surv(cclinic_data$time, cclinic_data$death) ~ cclinic_data$cluster)
colors = colors[groups]

pl = ggsurvplot(fit,
                censor = F,
                data = cclinic_data,
                risk.table = T, 
                risk.table.height = 0.3, 
                linetype = 1, 
                pval = TRUE,
                palette = colors,
                legend.labs = names(colors),
                ncensor.plot = FALSE)
print(pl)


# As we can see, even if the p-value is significant, survival curves cross, 
# suggesting that the assumption of constant hazards employed by Cox regression 
# is not valid in these patients. 
# 
# Thus, the observed p-value is not very interesting, and it is better to focus
# on a pairwise test. 
# 
# We can render this formal with a log-rank test that compares the survival
# curves for pairs of groups

# We write a function to repeat some operations
pairwise.test = function(clinical_data, groups) 
{
  clinical_data = clinical_data[which(clinical_data$cluster%in%groups), ]
  
  p = survdiff(
    Surv(clinical_data$time, clinical_data$death) ~ clinical_data$cluster, 
    data = clinical_data)
  
  # the p-valye is obtained as usual
  p.val = 1 - pchisq(p$chisq, length(p$n) - 1)
  
  return(p.val)
}

# We create a matrix of tests
test = matrix(0, length(groups),length(groups))
test = as.data.frame(test)

colnames(test) = sort(unique(cclinic_data$cluster))
rownames(test) = sort(unique(cclinic_data$cluster))

# ... and compute pairwise the values
for(i in rownames(test)){
  for (j in setdiff(colnames(test), i)) {
    test[i,j] = pairwise.test(cclinic_data, groups = c(i,j))  
  }
}

# These are the p-values
test

# These are the significant one at level 0.05
test < 0.05


# We conclude creating some plots of the most interesting findings.
corrplot::corrplot(as.matrix(test),
         col = colorRampPalette(RColorBrewer::brewer.pal(8, 'Greens'))(200),
         type = "upper", 
         p.mat = as.matrix(test), 
         number.font = 9,
         number.cex = 9,
         number.digits = 2,
         sig.level = -1,
         tl.cex = 1.5,
         insig = "p-value",
         tl.col = "black",
         cl.cex = 1,
         is.corr = F, 
         diag = F)

# Plots for separate groups
for (i in rownames(test)) {
  for (j in setdiff(colnames(test), i)) {
    if (test[i, j] < 0.05) 
    {
      cd1 = cclinic_data[which(cclinic_data$cluster%in%i), ] 
      cd2 = cclinic_data[which(cclinic_data$cluster%in%j), ] 
      
      cd = rbind(cd1,cd2)
      
      fit = survfit(Surv(cd$time, cd$death)~cd$cluster)
      # Adjust color label
      if(j>i){
      pl = ggsurvplot(fit,
                      data = cd,  
                      linetype=c(1,1),
                      pval = TRUE,
                      palette = colors[c(i,j)],
                      legend.labs = names(colors[c(i,j)]),
                      risk.table=F,
                      censor = T, 
                      conf.int = TRUE,
                      ncensor.plot = F)
      print(pl)
      }
    }
  }
}

