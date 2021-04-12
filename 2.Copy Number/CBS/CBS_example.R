require(dplyr)
require(ggplot2)

# 2 Poissons, mean 4 and 1
s1 = rpois(n = 100, lambda = 4) %>%
  as_tibble() %>%
  mutate(i = row_number())

s2 = rpois(n = 200, lambda = 1) %>%
  as_tibble() %>%
  mutate(i = row_number() + nrow(s1))

x = bind_rows(s1, s2)

x %>%
  ggplot(aes(x = i , y = value)) +
  geom_point() +
  geom_vline(xintercept = nrow(s1))

# The original changepoint detection strategy comes from a 1975 paper by Sen and
# Srivastava in the Annals of Statistics. If  D  is the array of data, of length  n,
# we consider the partial means  μ i  and  μ ′ i  of the first  i  and last  n − i
# elements of  D . We locate  i  such that  | μ i − μ ′ i |  is maximal and we apply
# a  t -test to determine if this difference is significant. If so, we mark this as a change point.

n = x %>% nrow

mean_cbs = function(i)
{
  mu_i = x %>%
    filter(row_number() < !!i) %>%
    pull(value) %>%
    mean
  
  mu_ip = x %>%
    filter(row_number() > n - i) %>%
    pull(value) %>%
    mean
  
  (mu_i - mu_ip) %>% abs
}

means = sapply(2:nrow(x), mean_cbs) %>%
  as_tibble()  %>%
  mutate(i = row_number())

ggplot(means, aes(x = i, y = value)) +
  geom_point() +
  geom_vline(xintercept = nrow(s1))

ttest_cbs = function(i)
{
  mu_i = x %>%
    filter(row_number() < !!i) %>%
    pull(value)
  
  mu_ip = x %>%
    filter(row_number() > n - i) %>%
    pull(value)
  
  if (length(mu_i) == 1)
    return(1)
  if (length(mu_ip) == 1)
    return(1)
  
  ttest = t.test(mu_i, mu_ip, alternative = 'greater')
  ttest$p.value
}

pvals = sapply(2:nrow(x), ttest_cbs) %>%
  as_tibble()  %>%
  mutate(i = row_number())

means_pvals = means %>% 
  rename(mean = value) %>% 
  full_join(pvals, by = 'i') %>%  
  rename(pval = value)

cut_p = (0.001)/nrow(means_pvals)

cutoff = means_pvals %>% 
  filter(pval < cut_p) %>% 
  arrange(desc(mean)) %>% 
  slice(1) %>% 
  pull(i)

ggplot(means_pvals, aes(
  x = i,
  y = mean,
  color = pval < cut_p)
  ) +
  geom_point() +
  geom_vline(xintercept = nrow(s1)) +
  geom_vline(xintercept = cutoff, color = 'red')

x %>%
  ggplot(aes(x = i , y = value)) +
  geom_point() +
  geom_vline(xintercept = nrow(s1)) +
  geom_vline(xintercept = cutoff, color = 'red')

# PSCBS
# 
# https://cran.r-project.org/web/packages/PSCBS/vignettes/CBS.pdf
install.packages('PSCBS')
require(PSCBS)

data <- PSCBS::exampleData("paired.chr01")
data <- data[, c("chromosome", "x", "CT")]
colnames(data)[3] <- "y"
data <- PSCBS::dropSegmentationOutliers(data)

gaps <- findLargeGaps(data, minLength = 1e+06)
knownSegments <- gapsToSegments(gaps)

fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
plotTracks(fit)

# Make x like data
data %>% head

# Run it on data, again
data = x %>% 
  mutate(chromosome = 1, x = i, y = value) %>% 
  dplyr::select(chromosome, x, y)
  
data <- PSCBS::dropSegmentationOutliers(data)

gaps <- findLargeGaps(data, minLength = 1e+06)
knownSegments <- gapsToSegments(gaps)

fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
plotTracks(fit)

# Add third Poisson
s3 = rpois(n = 120, lambda = 2) %>%
  as_tibble() %>%
  mutate(i = row_number() + nrow(x))

x = x %>% bind_rows(s3)

# Run again
data = x %>% 
  mutate(chromosome = 1, x = i, y = value) %>% 
  dplyr::select(chromosome, x, y)

data <- PSCBS::dropSegmentationOutliers(data)

gaps <- findLargeGaps(data, minLength = 1e+06)
knownSegments <- gapsToSegments(gaps)

fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
plotTracks(fit)

