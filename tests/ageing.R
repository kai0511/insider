####################################################################################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
####################################################################################################################

require(insider)

# options(error = dump_and_quit())

num_factors <- 23

# regularization for iMF
lambda <- 10
alpha <- 0.4

setwd("../results/ageing")
load('~/data/multi_dimensional_datasets/ageing_dataset_annotated_with_phenotypes_filtered.RData')

dataset[is.na(dataset)] <- 0
dataset <- dataset[,-1]

end_idx <- 3
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)
confounders <- as.matrix(dataset[ ,1:end_idx])

object <- insider(data, as.matrix(confounders), as.integer(c(1,2)), global_tol = 1e-10, sub_tol = 1e-5, tuning_iter = 30)
# object <- tune(object, as.integer(num_factors), lambda = seq(1, 20, 1), alpha = c(0.2, 0.3, 0.4, 0.5))
object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = paste0("insider_ageing_R", num_factors, "_fitted_object.RData"))
