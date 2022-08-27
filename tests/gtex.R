#######################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
#######################################################
require(insider)

num_factors <- 13

# regularization for iMF
lambda <- 44.56  # original value 33.7777777777778
alpha <- 0.4 # original 0.35


setwd('../results/gtex/')

load("~/data/multi_dimensional_datasets/gtex_brain_sampled_expression_including_phenotype.RData")

# dataset <- as.matrix(dataset)
# dataset[is.na(dataset)] <- 0
data <- log2(as.matrix(dataset[,-c(1,2,3)]) + 1)
confounders <- as.matrix(dataset[,2:3])
colnames(confounders) <- c('gender', 'structure')


object <- insider(data, as.matrix(confounders), global_tol = 1e-10)
object <- tune(object, latent_dimension = as.integer(seq(10, 20, by = 1)), lambda = 1, alpha = 0.1) 
# object <- fit(object, latent_dimension = as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = "insider_brainspan_fitted_object.RData")
