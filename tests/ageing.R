####################################################################################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
####################################################################################################################

require(insider)

# options(error = dump_and_quit())

num_factors <- 26

# regularization for iMF
lambda <- 15
alpha <- 0.2

setwd("../results/ageing")
load('~/data/multi_dimensional_datasets/ageing_dataset_annotated_with_phenotypes_filtered.RData')

dataset[is.na(dataset)] <- 0
dataset <- dataset[,-1]

end_idx <- 3
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)
confounders <- as.matrix(dataset[ ,1:end_idx])

object <- insider(data, as.matrix(confounders), global_tol = 1e-10, tuning_iter = 30)
# object <- tune(object, as.integer(num_factors), lambda = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50), alpha = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = "insider_ageing_fitted_object.RData")
