#install the simulator
install.packages("Rcpp")   # C++ compiler for r
library(Rcpp)

devtools::install_github("T-Heide/TEMULATOR")   # T.Heide simulator

library(TEMULATOR) # package containing the simulator and datasets
library(ggplot2)
library(dplyr)

#We look at some VAF spectrum used in the Mobster paper

data("mobster_simulations", package="TEMULATOR")

# example of Non-neutral simulation

# get info of the simulations
print(mobster_simulations[[1]]) 

# plot the VAF spectrum
plot(mobster_simulations[[1]])     

# Selection of suitable simulations

data("mobster_summary_of_simulations", package="TEMULATOR")

mobster_summary_of_simulations$group =
  case_when(mobster_summary_of_simulations$subclone_fraction < 0.05 ~ "neutral_low_frac",
            mobster_summary_of_simulations$subclone_fraction > 0.9 ~ "neutral_high_frac",
            between(mobster_summary_of_simulations$subclone_fraction, 0.2, 0.8)  ~ "detectable")

mobster_summary_of_simulations_subset =
  mobster_summary_of_simulations %>%
  filter(n_mutations_before_insertions >= 50) %>%
  filter(!is.na(group)) %>%
  filter(!is.na(assigned_id)) %>%
  arrange(subclone_fraction)


ids_per_group =
  split(mobster_summary_of_simulations_subset$assigned_id,
        mobster_summary_of_simulations_subset$group)

sapply(ids_per_group, length)


# Example 1: Neutral evolution with low subclone fraction

set_1= "neutral_low_frac"
low_fraction_simultations=mobster_simulations[ids_per_group[[set_1]]] # list of TEMULATOR object
plot(low_fraction_simultations$e6a31b125cd60319a05deced35775d02651a4e3c) # choose a simulation


# Example 2: Neutral evolution with high subclone fraction

set_2 = "neutral_high_frac"
high_fraction_simultations= mobster_simulations[ids_per_group[[set_2]]] # list of TEMULATOR object
plot(high_fraction_simultations$f7bc3f781663dd18d5e08ab3df198a65714a708c) # choose a simulation

# Example 3: VAF spectrum with detectable subclone

set_3 = "detectable"
detectable_simulations = mobster_simulations[ids_per_group[[set_3]]] # list of TEMULATOR object
plot(detectable_simulations$`73df1a4e7d8afe464aa218e276f3fac6fa40ded6`)  # choose a simulation


# Recreate a simulation with TEMULATOR

# Simulation parameters:
birthrates=c(1.0, 2.0)
mutation_rate = 16   # mutation rate per cell doubling (m)
death_rate = 0.2     # death rate (μ)
end_time = 782830 # number of reactions until the simulation in finished
sc_st =  1000    # subclone start time

# Sequencing parameters:
depth_model = 2   # over-dispersed beta-binomial distribution, dispersion parameter ρ=0.0
n_clonal = 500    # number of clonal mutations (N_clonal)
depth = 120       # 120x simulated coverage (C̅)
min_vaf = 0.05    # minimum VAF to accept a variant

simulation_result =
  simulateTumour(
    birthrates = birthrates,
    deathrates = rep(death_rate, 2),
    mutation_rates = rep(mutation_rate, 2),
    clone_start_times = c(0, sc_st),
    fathers = c(0, 0),                  # CLONAL -> SC1 -> SC2, linear evolution 
    simulation_end_time = end_time,
    seed = seed,
    depth_model = depth_model,
    depth = depth,
    min_vaf = min_vaf,
    number_clonal_mutations = n_clonal,
    verbose=TRUE
  )

plot(simulation_result)

# more at https://t-heide.github.io/TEMULATOR/index.html

