##############################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
##############################################

# load('~/data/zscore_matrix_selected.RData')
# 
# dataset[is.na(dataset)] <- 0
# splitted <- ratio_splitter(dataset)
# trainset <- splitted$trainset
# testset <- splitted$testset
# train_indicator <- (trainset != 0)
# test_indicator <- (testset != 0)
# 
# names <- rownames(trainset)
# grouping <- t(unname(sapply(names, function(s) split_str(s))))
# # grouping <- grouping[grouping[ , 1] == single_disease, ]
# confounders <- as.data.frame(grouping, stringsAsFactors = T)
# colnames(confounders) <- c('disease', 'tissue')
# num_disease <- nlevels(confounders[[1]])
# num_tissue <- nlevels(confounders[[2]])
# levels(confounders[[1]]) <- seq(num_disease)
# levels(confounders[[2]]) <- seq(num_tissue)
# confounders$disease <- as.integer(confounders[[1]])
# confounders$tissue <- as.integer(confounders[[2]])
# 
# num_factors <- 50
