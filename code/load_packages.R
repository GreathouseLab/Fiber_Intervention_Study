# ============================================= #
# script: load_packages.R
# Project: Fiber Intervention Study
# Author(s): L. Greathouse et al.
# ============================================= #
# Date Created: 2019-12-10
# Date Modified: 2019-12-10
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading all necessary
#   R packages
#
# No output - just loading packages into the
#   environment
# ============================================= #
# Set up directory and libraries
rm(list=ls())
# list of packages
packages <- c("phyloseq",
              "tidyverse", "readr", "readxl",
              "data.table", "dplyr","ggplot2",
              "kableExtra", "xtable")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Load packages
lapply(packages, library, character.only = TRUE)

w.d <- getwd()

