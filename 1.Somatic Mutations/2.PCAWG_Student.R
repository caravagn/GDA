# AT HOME
# 
# Download the PCAWG cohort data
#
# -> https://dcc.icgc.org/releases/PCAWG/
# -> https://dcc.icgc.org/releases/PCAWG/consensus_snv_indel
# -> https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz
#
# ~2700 WGS pan-cancer, ~1GB of somatic mutations calls.
# 
# Write your own parsers for the data, load it and store it in some convenient way (eg, an rds per sample). To complete
# this you will have to complete a few more lectures, but you can start pooling mutation data after this lecture

##########################################################################################
# Lecture - Somatic Mutation Calling
##########################################################################################

# I assembled some PCAWG cases for you (somatic calls)
PCAWG_example = "https://www.dropbox.com/s/stkkocdjymmtcp8/PCAWG_somatic_mutations.zip?raw=1"

# Download, load and cancel data
download.file(PCAWG_example, "PCAWG_somatic_mutations.zip")

unzip("PCAWG_somatic_mutations.zip")

file.remove("PCAWG_somatic_mutations.zip")

########################################################################################## 
# Download, unzip and cancel data, then load all the data in a tibble. 
# There should be ~500.000 mutations, add tumor type joining data file Types.csv, 
# and add a categorical flag "SNV", "indel" based on the type of mutations identified 
# (use reference and alternative alleles to determine that)
##########################################################################################
all_files = list.files('./samples/', full.names = TRUE)

##########################################################################################
# Visualise the distribution of SNVs/indels across tumour types, with a piechart/ barplot
##########################################################################################

##########################################################################################
# Compute mutational burden per sample, split by type of mutation
##########################################################################################

##########################################################################################
# Define the substitutions for an SNV as C>T (REF>ALT), compute their distribution per cancer type.
# Which tumour has more mutations? Can you guess why?
##########################################################################################

##########################################################################################
# Compute the mutational burden per tumour type, sort types by that, and use a scatter
# where you can annotate each sample by median depth of sequencing
##########################################################################################

##########################################################################################
# Compute which genes has the highest number of "Missense_Mutation" flags. Do you
# see any gene you ever heard of?
##########################################################################################

##########################################################################################
# Reproduce this plot: 
# 
# https://github.com/cBioPortal/icebox/issues/78 
##########################################################################################

##########################################################################################
# Remake these plots with all the PCAWG cohort. Note that if you download the data you
# will have to process data columns etc (I polished those for you now).
##########################################################################################

