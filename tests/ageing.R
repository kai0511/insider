####################################################################################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
####################################################################################################################

require(insider)

options(error = dump_and_quit())

num_factors <- 23

# regularization for iMF
lambda <- 7
alpha <- 0.2

# setwd("~/data/multi_dimensional_datasets/result/ageing/")
load('~/data/multi_dimensional_datasets/ageing_dataset_annotated_with_phenotypes_filtered.RData')
dataset[is.na(dataset)] <- 0
dataset <- dataset[,-1]
end_idx <- 3

confounders <- dataset[ ,1:end_idx]
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)



