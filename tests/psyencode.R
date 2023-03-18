##############################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
##############################################
require(insider)

# options(error = dump_and_quit)

num_factors <- 3

# regularization for insider
lambda <- 120  # original value 33.7777777777778
alpha <- 0.9 # original 0.35

setwd('~/git_repositories/insider/results/PsyEncode/')
load('~/data/PsyEncode/PEC_Gene_expression_matrix_for_insider.RData')
dataset[is.na(dataset)] <- 0

end_idx <- 3
confounders <- dataset[ ,1:end_idx]
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)

object <- insider(data, as.matrix(confounders), global_tol = 1e-10)
# object <- tune(object, latent_dimension = as.integer(num_factors), lambda = seq(90, 500, by = 10), alpha = 0.9)
object <- fit(object, latent_dimension = as.integer(num_factors), lambda = lambda, alpha = alpha, partition = 0)
save(object, file = paste0("insider_psyencode_R", num_factors, "_fitted_object.RData"))
 
