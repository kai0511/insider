##############################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
##############################################
<<<<<<< HEAD
# options(error = dump_and_quit())
=======
options(error = dump_and_quit)
>>>>>>> c3b3f6b3bf67ec0b08cdf83fb4df61f5bd8ef18a

num_factors <- 32 

# regularization for iMF
lambda <- 12.22  # original value 33.7777777777778
alpha <- 0.6 # original 0.35

# load('~/data/zscore_matrix_selected.RData')
load('~/data/multi_dimensional_datasets/brainspan_dataset_annotated_fitered.RData')
## load('~/data/multi_dimensional_datasets/aging_dataset_annotated_with_phenotypes_filtered.RData')
dataset[is.na(dataset)] <- 0
# dataset <- dataset[,-1]
end_idx <- 2

confounders <- dataset[ ,1:end_idx]
confounders[, 1] <- confounders[, 1] - 1
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)
 
