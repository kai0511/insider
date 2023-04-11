##############################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
##############################################
require(insider)

# options(error = dump_and_quit)

num_factors <- 23 

# regularization for iMF
lambda <- 6  # original value 33.7777777777778
alpha <- 0.4 # original 0.35

setwd('../results/brainspan/')
load('~/data/multi_dimensional_datasets/brainspan_dataset_annotated_fitered.RData')
dataset[is.na(dataset)] <- 0

end_idx <- 2
confounders <- dataset[ ,1:end_idx]
confounders[, 1] <- confounders[, 1] - 1
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)

object <- insider(data, as.matrix(confounders), global_tol = 1e-10)
# object <- tune(object, latent_dimension = as.integer(num_factors), lambda = seq(1, 50, by = 5), alpha = seq(0.1, 0.6, by = 0.1))
object <- fit(object, latent_dimension = as.integer(num_factors), lambda = lambda, alpha = alpha, partition = 0)
# save(object, file = paste0("insider_brainspan_R", num_factors, "_fitted_object.RData"))
 
