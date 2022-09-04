#######################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
#######################################################
require(insider)

num_factors <- 12

# regularization for iMF
lambda <- 60  # original value 33.7777777777778
alpha <- 0.5 # original 0.35


setwd('../results/gtex/')

load("~/data/multi_dimensional_datasets/gtex_brain_sampled_expression_including_phenotype.RData")

# dataset <- as.matrix(dataset)
# dataset[is.na(dataset)] <- 0
data <- log2(as.matrix(dataset[,-c(1,2,3)]) + 1)
confounders <- as.matrix(dataset[,2:3])
colnames(confounders) <- c('gender', 'structure')


object <- insider(data, as.matrix(confounders), as.integer(c(1,2)), global_tol = 1e-10)
# object <- tune(object, latent_dimension = as.integer(12), lambda = seq(30, 60, by = 5), alpha = c(0.3, 0.4, 0.5, 0.6)) 
object <- fit(object, latent_dimension = as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = "insider_gtex_fitted_object.RData")
