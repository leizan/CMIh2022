#This function load all the packages and functions needed for the discretization.

# required for algorithm
require("parallel")
require("SCCI")
require(dplyr)
library("MASS")
library(foreach)
library(doParallel)
library(extraDistr)
library(parallelDist)
library(qlcMatrix )

# required for experiments
require("pcalg")


Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")

source("MS-R.R")
source("experiments/syntheticGenerator.R")
source("mixedCmiIEstimator.R")
source("experiments/data.R")
source("method.R")


source("algorithm/generate_candidate_cuts.R")
source("algorithm/refine_candidate_cuts_exclude_empty_bins.R")
source("algorithm/modified_log.R")
source("algorithm/iterative_cmi_greedy_flexible_parallel.R")

source("algorithm/multi_hist_splitting_seed_based_simple.R")
source("algorithm/oned_hist_iterative.R")
source("algorithm/CMI_estimates.R")
source("algorithm/CMI_pvals.R")
source("algorithm/utils.R")
