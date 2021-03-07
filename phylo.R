require(dplyr)
require(ggplot2)

# Made up sequence distance computation
DNA_alphabet = c("A", "C", "T", "G")

sequence_length = 20

sequences = 15

# Samples
samples = lapply(1:sequences, function(x)
  sample(
    DNA_alphabet,
    sequence_length,
    prob = c(.8, .1, .1, .2),
    replace = TRUE
  ))
names(samples) = paste0("S", 1:sequences)

# Sample-to-sample distance (indicator function)
f_dist = function(samples, i, j) {
  (samples[[i]] != samples[[j]]) %>% sum
}

# Matrix of pairwise distances
M_dist = function(samples) {
  n = samples %>% length
  M = matrix(0, nrow = n, ncol = n)
  colnames(M) = rownames(M) = names(samples)
  
  for (i in colnames(M))
    for (j in colnames(M))
      M[i, j] = f_dist(samples, i, j)
  
  M = M / max(M)
  
  M %>%
    as.dist()
}

M_dist(samples)

my_nj <- ape::nj(M_dist(samples))
plot(my_nj, "unrooted")
plot(my_nj)

# ggtree plot (alternative plotting)
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("ggtree")
#
# library("ggtree")
#
# ggtree(my_nj, layout="equal_angle", size=0.5, linetype=1)

# Adenocarcinoma Set06
multi_region = evoverse.datasets::MSEQ_CRC_ADENOCARCINOMA_SET6

multi_region$mutations = multi_region$mutations %>%
  filter(region == 'exonic')

n = multi_region$samples %>% length
M = matrix(0, nrow = n, ncol = n)
colnames(M) = rownames(M) = multi_region$samples

for (i in colnames(M))
  for (j in colnames(M))
  {
    s1 = multi_region$mutations[[paste0(i, '.VAF')]] > 0
    s2 = multi_region$mutations[[paste0(j, '.VAF')]] > 0
    
    M[i, j] = (s1 != s2) %>% sum
  }

M = M / max(M)

M = M %>%
  as.dist()

my_nj <- ape::nj(M)
plot(my_nj, "unrooted")
plot(my_nj)

ggtree(my_nj,
       color = "indianred3",
       size = .5,
       layout = "roundrect") +
  geom_tiplab(color = 'purple4') +
  geom_nodepoint(color = "#b5e521",
                 alpha = 1 / 4,
                 size = 10) +
  geom_tippoint(color = "steelblue",
                shape = 8,
                size = 3) +
  theme_tree("gainsboro")


ggtree(my_nj,
       color = "indianred3",
       branch.length = 'none',
       layout = 'circular') +
  geom_tiplab(color = 'purple4') +
  geom_tippoint(color = "steelblue",
                shape = 8,
                size = 3) +
  geom_nodepoint(color = "#b5e521",
                 alpha = 1 / 4,
                 size = 10) +
  theme_tree("gainsboro")


# Inflating fake selection
M = M %>% as.matrix()
M[2:nrow(M), 1] = M[2:nrow(M), 1] + 3

my_nj <- ape::nj(M)

plot(my_nj, "unrooted")
plot(my_nj)




require(dplyr)

# Plot mutliple trees with aligned tips
trees <-
  list(read.tree(text = "((a:1,b:1):1.5,c:2.5);"),
       read.tree(text = "((a:1,c:1):1,b:2);"))

ggdensitree(trees) + geom_tiplab()

# Plot multiple trees with aligmned tips with tip labls and separate tree colors
trees.fort <-
  list(trees[[1]] %>% fortify %>% mutate(tree = "a"),
       trees[[2]] %>% fortify %>% mutate(tree = "b"))

ggdensitree(trees.fort, aes(colour = tree)) + geom_tiplab(colour = 'black')


# Read example data
trees <-
  read.tree(system.file("examples", "ggdensitree_example.tree", package =
                          "ggtree"))

# Compute OTU
grp <-
  list(
    A = c("a.t1", "a.t2", "a.t3", "a.t4"),
    B = c("b.t1", "b.t2", "b.t3", "b.t4"),
    C = c("c.t1", "c.t2", "c.t3", "c.t4")
  )
trees <- lapply(trees, groupOTU, grp)

# Plot multiple trees colored by OTU
ggdensitree(trees, aes(colour = group), alpha = 1 / 6) + scale_colour_viridis_d(option = 'D')


myBoots <-
  boot.phylo(my_nj, M %>% as.matrix, function(xx)
    nj(M), B = 10,  mc.cores = 6, trees = TRUE)
ggdensitree(myBoots$trees, alpha = 1 / 6) + scale_colour_viridis_d(option = 'D')

# TRACERx
download.file(
  "https://github.com/caravagn/revolver.misc/blob/master/vignette_TRACERx_Hanjani_et_al/%5BData%5D%20TRACERx.appendix2.Rdata?raw=true",
  destfile = 'Tracerx.RData'
)

load("Tracerx.RData")

data = data %>% as_tibble()

patient = "CRUK0024"

get_pat = function(patient)
{
  patient_data = data %>% filter(SampleID == patient,!is.na(ObsPyCloneCCF))
  
  df = NULL
  
  for (i in 1:nrow(patient_data))
  {
    w = revolver::CCF_parser(patient_data$ObsPyCloneCCF[i]) %>% as.numeric
    w[w <= .2] = 0 %>% as.integer()
    w[w > .2] = 1 %>% as.integer()
    df = rbind(df, w)
  }
  
  df = df %>% as_tibble() %>% distinct
  colnames(df) = paste0("R", 1:ncol(df))
  
  df= df[complete.cases(df), ]
  df %>% as.matrix()
}


patient_data = get_pat(patient)

distance_M = function(x)
{
  n = x %>% ncol()
  M = matrix(0, nrow = n, ncol = n)
  colnames(M) = rownames(M) = x %>% colnames
  
  for (i in colnames(M))
    for (j in colnames(M))
    {
      s1 = x[, i] > 0
      s2 = x[, j] > 0
      
      M[i, j] = (s1 != s2) %>% sum
    }
  
  M = M / max(M)
  
  M %>%
    as.dist()
}


my_nj <- ape::nj(distance_M(patient_data))
plot(my_nj, "unrooted")
plot(my_nj)

patient_data_matrix = patient_data %>% as.matrix() %>% t

distance_M(x = patient_data_matrix %>% t)

myBoots <- boot.phylo(my_nj,
                      patient_data_matrix,
                      function(xx) {
                        xxrd = xx %>% t
                        xxrd = xxrd[, sort(xxrd %>% colnames)]
                        
                        print(xxrd)
                        nj(distance_M(xxrd))
                      },
                      B = 30,
                      mc.cores = 6,
                      trees = TRUE)$trees 

ggdensitree(myBoots, alpha = 1 / 6) + scale_colour_viridis_d(option = 'D')

ggdensitree(myBoots$trees[1:2], 
            alpha=.3,   
            layout = "slanted",
            colour='steelblue') + 
  geom_tiplab(size=3) + 
  xlim(0, 2)


# Plot mutliple trees with aligned tips
trees <- list(read.tree(text="((a:1,b:1):1.5,c:2.5);"), read.tree(text="((a:1,c:1):1,b:2);"));
ggdensitree(trees) + geom_tiplab()

trees.fort <- list(trees[[1]] %>% fortify %>% mutate(tree="a"), trees[[2]] %>% fortify %>% mutate(tree="b"));
ggdensitree(myBoots_trees, aes(colour=tree)) + geom_tiplab(colour='black')


plot(myBoots_trees[[1]])


# Generate example data
set.seed(1)
trees <- rmtree(5, 10)
time.trees <- lapply(1:length(trees), function(i) {
  tree <- trees[[i]]
  tree$tip.label <- paste0("t", 1:10)
  dates <- estimate.dates(tree, 1:10, mu=1, nsteps=1)
  tree$edge.length <- dates[tree$edge[, 2]] - dates[tree$edge[, 1]]
  fortify(tree) %>% mutate(tree=factor(i, levels=as.character(1:10)))
})

# Plot multiple trees with aligned tips from muliple time points
ggdensitree(myBoots_trees, aes(colour=tree), tip.order=paste0("t", 1:10)) + geom_tiplab(colour='black')


