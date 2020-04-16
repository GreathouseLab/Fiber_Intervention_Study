# ============================================= #
# script: get_data.R
# Project: Fiber Microbiome Study
# Author(s): R.N. Padgett & L. Greathouse
# ============================================= #
# Data Created: 2019-10-16
# Date Modified: 2020-04-10
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading and formating
#   the data file for use in analyses.
#
# ============================================= #

# need to make sure the microbiome code is read in
source("code/microbiome_statistics_and_functions.R")

# read in analysis
# from data-for-analaysis tab of
analysis_data <- read_xlsx("data/analysis-data/Data_Fiber_2020_02_06.xlsx", sheet="DataforAnalysis",na = ".")

# supplement IDS
# *highest intake = 7, supplement
# highest cookie intake = 3
# Gut health = 1 = no issues
# 2 = some discomfort
# 3 = sig discomfort/constipation
# 4 = sig discomfort/diarrhea


# get microbiome data
biom_file  <- import_biom("data/analysis-data/OTU_Table.biom")
tree_file <- read_tree("data/analysis-data/OTU_Table.tre")
meta_data <- read.table("data/analysis-data/Metadata file_nomiss_microbiomeIDs.txt", sep="\t", header=T)
  colnames(meta_data)[1] <- "ID"
  rownames(meta_data) <- meta_data$ID
  meta_data$SubjectID2 <- as.numeric(as.factor(meta_data$SubjectID))
  meta_data$Week <- as.factor(meta_data$Week)
  ## merge meta_dat with analysis data
  meta_data <- merge(meta_data, analysis_data, all=T)
meta <- sample_data(meta_data)
sample_names(meta) <- meta_data$ID
phylo_data0 <- merge_phyloseq(biom_file, tree_file, meta)
colnames(phylo_data0@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# phylo_data <- rarefy_even_depth(phylo_data, sample.size = 7004, rngseed=20191210)
# phylo object ready for phyloseq related analyses (alphe, beta, etc..)

# now, extract the information from .biom/phyloseq for other analyses
otus <- psmelt(biom_file) # uses the .biom file
otus <- otus[ otus$Sample %in% meta_data$ID, ]
otus <- reshape(otus, idvar = "OTU", timevar = "Sample", direction = "wide")
OTU <- otus$OTU
Kingdom <- otus[, paste0('Rank1.', meta_data$ID[1])]
Phylum <- otus[, paste0('Rank2.', meta_data$ID[1])]
Class <- otus[, paste0('Rank3.', meta_data$ID[1])]
Order <- otus[, paste0('Rank4.', meta_data$ID[1])]
Family <- otus[, paste0('Rank5.', meta_data$ID[1])]
Genus <- otus[, paste0('Rank6.', meta_data$ID[1])]
# Next get just the observed counts of the OTUs
observed_counts <- otus[, grep("Abundance",colnames(otus)) ]
i <- 1
for(i in 1:length(colnames(observed_counts))){
  colnames(observed_counts)[i] <- substring(colnames(observed_counts)[i], 11)
}
# OTU Names object
otu.name <- cbind(Kingdom,Phylum,Class,Order,Family,Genus)
rownames(otu.name) <- OTU
# Abundance List Objects
abund.list <- list( cbind(Kingdom,observed_counts),
                    cbind(Phylum,observed_counts),
                    cbind(Class,observed_counts),
                    cbind(Order,observed_counts),
                    cbind(Family,observed_counts),
                    cbind(Genus,observed_counts),
                    cbind(OTU,observed_counts))
names(abund.list) <- c("Kingdom","Phylum","Class","Order","Family","Genus","OTU")
i <- 1
for(i in 1:7){
  abund.list[[i]] <- abundance_list_create(abund.list[[i]],abund.list[[i]][,1])
}


microbiome_data <- list(otu.tab = abund.list$OTU,
                        otu.name,abund.list,
                        meta.dat = meta_data,
                        tree = tree_file)
names(microbiome_data) <- c("otu.tab", "otu.name","abund.list",
                                  "meta.dat","tree")

# remove unnecessary items
remove(abund.list, analysis_data, biom_file, meta, meta_data, observed_counts, otu.name, otus, tree_file, Class, Family, Genus, i, Kingdom, new.packages, Order, OTU,packages, Phylum)
