##############################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
##############################################
require(insider)

# options(error = dump_and_quit)

num_factors <- 19 

# regularization for iMF
lambda <- 16  # original value 33.7777777777778
alpha <- 0.5 # original 0.35

setwd('../results/brainspan/')
load('~/data/multi_dimensional_datasets/brainspan_dataset_annotated_fitered.RData')
dataset[is.na(dataset)] <- 0

end_idx <- 2
confounders <- dataset[ ,1:end_idx]
confounders[, 1] <- confounders[, 1] - 1
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)

object <- insider(data, as.matrix(confounders))
object <- tune(object, latent_dimension = as.integer(19), lambda = c(10, 12, 14, 16, 18, 20, 22, 24, 26, 30), alpha = c(0.4, 0.5, 0.6, 0.7))
# object <- fit(object, latent_dimension = as.integer(c(17, 18, 19, 20)), lambda = lambda, alpha = alpha)
 
