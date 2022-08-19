####################################################################################################################
# This file models 3-D matrix factorization with additive matrix factorization using L1 penalty.
# Alternating Minimization and Coordinate Descent are employed to optimize this problem.
# For details see my note, derivartion for matrix factorization.
# This version is a mature alpha version, and some improvement is in progess.
####################################################################################################################

require(dplyr)
require(methods)
require(parallel)
require(Rcpp)

# load global configuration
source('./global_conf.R')
source('./imputation_functions.R')
# source('./optimization_functions.R')
sourceCpp("./cd_screening.cpp")

# dump objects for easy debugging
dump_and_quit <- function() {
  # Save debugging info to file last.dump.rda
  dump.frames(to.file = TRUE)
  # Quit R with error status
  q(status = 1)
}
options(error = dump_and_quit)

predict <- function(disease_factor, tissue_factor, donor_factor, column_factor, cfd1, cfd2, cfd3) {
    # disease_matrix <- disease_factor %*% column_factor
    # tissue_matrix <- tissue_factor %*% column_factor
    # donor_matrix <- donor_factor %*% column_factor

    predictions <- (disease_factor[cfd1, ] + tissue_factor[cfd2, ] + donor_factor[cfd3, ] + disease_factor[cfd1, ] * tissue_factor[cfd2, ]) %*% column_factor

    return(predictions)
}

compute_loss <- function(disease_factor, tissue_factor, donor_factor, column_factor, sum_residual, verbose = T){

    # l2 penalty 
    row_reg_loss <- lambda * (norm(tissue_factor, "F")^2 + norm(disease_factor, "F")^2 + norm(donor_factor, "F")^2)
    col_reg_loss <- lambda * (1 - alpha) * norm(column_factor, "F")^2
    
    # l1 penalty
    l1_reg_loss <- alpha * lambda * sum(abs(column_factor))

    if(verbose){
        cat(paste0('total_residual:\t', sum_residual/2, 
                   '; \nrow_reg_loss:\t', row_reg_loss/2, 
                   '; \ncol_reg_loss:\t', col_reg_loss/2, 
                   '; \nl1_reg_loss:\t', l1_reg_loss, "\n"))
    }
    loss <- sum_residual/2 + row_reg_loss/2 + col_reg_loss/2 + l1_reg_loss
    return(loss)
}

evaluate_model <- function(trainset, testset, prediction, iter = NA, verbose = T, debug = F){

    train_diff <- (trainset - prediction)
    sum_residual <- sum((train_diff[train_indicator])^2)
    train_rmse <- sqrt(sum_residual/sum(train_indicator))

    test_diff <- ((testset - prediction)[test_indicator])^2
    test_rmse <- sqrt(mean(test_diff))

    if (verbose){
        train_info <- ifelse(!is.na(iter), 
                             paste0("ibiasedMF with L1_penalty iter ", iter,  ": train rmse = ", train_rmse), 
                             paste0("train rmse of ibiasedMF model: ", train_rmse))

        test_info <- ifelse(!is.na(iter), 
                            paste0("ibiasedMF with L1_penalty iter ", iter,  ": test rmse = ", test_rmse), 
                            paste0("test rmse of ibiasedMF model: ", test_rmse))
        cat(train_info, '\n')
        cat(test_info, '\n')
    }
    return(list(test_rmse = test_rmse, sum_residual = sum_residual))
}

save_model <- function(fitted_obj, iter){
    save(fitted_obj, file = paste0('Err_gene_expression_iMF_L1_penalty_',  iter, '.RData'))
}

train_amf_model <- function(trainset, confounders, num_factors = 50, num_iterations = 80000, verbose = T){
    # for disease tissue specific, confounders should have two columns, one for disease and another for tissue
    decay <- 1
    num_col <- ncol(trainset)
    col_loop <- seq(num_col)
    zero_idx <- (trainset == 0)

    cfd1 <- as.integer(confounders[[1]])	
    cfd2 <- as.integer(confounders[[2]])
    cfd3 <- as.integer(confounders[[3]])
    num_disease <- length(unique(cfd1))
    num_tissue <- length(unique(cfd2))
    num_donor <- length(unique(cfd3))

    
    disease_factor <- matrix(generate_rand_vector(num_factors * num_disease), ncol = num_factors)
    tissue_factor <- matrix(generate_rand_vector(num_factors * num_tissue), ncol = num_factors)
    donor_factor <- matrix(generate_rand_vector(num_factors * num_donor), ncol = num_factors)
    column_factor <- matrix(generate_rand_vector(num_factors * ncol(trainset)), nrow = num_factors)
    

    load("gene_expression_iMF_L1_penalty_23v1.RData") 
    disease_factor <- fitted_obj$disease_factor
    tissue_factor <- fitted_obj$tissue_factor
    donor_factor <- fitted_obj$donor_factor
    column_factor <- fitted_obj$column_factor

    indicator <- apply(train_indicator, 2, as.numeric)
    predictions <- predict(disease_factor, tissue_factor, donor_factor, column_factor, cfd1, cfd2, cfd3)
    metrics <- evaluate_model(trainset, testset, predictions)
    loss <- compute_loss(disease_factor, tissue_factor, donor_factor, column_factor, metrics$sum_residual)

    for(iter in seq(num_iterations)){
        last_loss <- loss
        if((iter %% 10 == 0) && verbose){
            cat("Iteration ", iter, '------------------------------------------\n')	
            cat('sum of squared error in iteration ', iter, ' before optimizating stage_factor: ', last_loss, '\n')
        }
        
        # fix all other parameters and update disease factor
        donor_matrix <- donor_factor %*% column_factor
        tissue_matrix <- tissue_factor %*% column_factor
        residual <- trainset - tissue_matrix[cfd2, ] - donor_matrix[cfd3, ]
        # residual[zero_idx] <- 0
        disease_factor <- optimize_row_l2(residual, indicator, disease_factor, tissue_factor, column_factor, cfd1, cfd2, lambda, n_cores = 16)
	    
        # fix all other parameters and update tissue factors
        disease_matrix <- disease_factor %*% column_factor
        residual <- trainset - disease_matrix[cfd1, ] - donor_matrix[cfd3, ]
        # residual[zero_idx] <- 0
        tissue_factor <- optimize_row_l2(residual, indicator, tissue_factor, disease_factor, column_factor, cfd2, cfd1, lambda, n_cores = 16)

        # fix all other parameters and update donor factors
        residual <- trainset - (disease_factor[cfd1, ] + tissue_factor[cfd2, ] + disease_factor[cfd1, ] * tissue_factor[cfd2, ])%*% column_factor
        # residual[zero_idx] <- 0
        # donor_factor <- optimize_row_l1(residual, stage_factor, tissue_factor, column_factor, cfd1, lambda, alpha, n_cores = 20, tol = 1e-12)
        donor_factor <- optimize_row_l2(residual, indicator, donor_factor, tissue_factor, column_factor, cfd3, cfd2, lambda, n_cores = 20, interaction = 0)

        # fix all other parameters and update column factor term
        row_factor <- disease_factor[cfd1, ] + tissue_factor[cfd2, ] + donor_factor[cfd3, ] + disease_factor[cfd1, ] * tissue_factor[cfd2, ]
        column_factor <- optimize_col(trainset, indicator, row_factor, column_factor, lambda, alpha, n_cores = 30, tol = 1e-7 * decay)

        if(iter %% 10 == 0){
            predictions <- predict(disease_factor, tissue_factor, donor_factor, column_factor, cfd1, cfd2, cfd3)
            metrics <- evaluate_model(trainset, testset, predictions)
            loss <- compute_loss(disease_factor, tissue_factor, donor_factor, column_factor, metrics$sum_residual)

            rmse <- metrics[['test_rmse']]
            cat(paste0('sum of squared error in iteration ', iter, ' after optimizing factor: ', rmse, '\n'))
            
            decay <- ifelse((last_loss - loss)/1000 >= 1, 1, signif((last_loss - loss)/1000, digits=1))
            if(is_converged(loss, last_loss, iter, learner = 'ibiasedMF_L1_penalty', thres = 1e-10)) {
                break
            }
            
            optimal_rmse <- rmse
        }
    }

    fitted_obj <- list (iter = iter,
                        trainset = trainset,
                        disease_factor = disease_factor,
                        tissue_factor = tissue_factor,
                        donor_factor = donor_factor, 
                        column_factor = column_factor, 
                        optimal_rmse = optimal_rmse)

    save(fitted_obj, file = paste0('gene_expression_iMF_L1_penalty_',  num_factors, 'v2.RData'))
    return(fitted_obj)
}

# regularization for iMF
lambda <- 7
alpha <- 0.2

setwd("~/data/multi_dimensional_datasets/result/ageing/")
load('~/data/multi_dimensional_datasets/ageing_dataset_annotated_with_phenotypes_filtered.RData')
dataset[is.na(dataset)] <- 0
dataset <- dataset[,-1]
end_idx <- 3

confounders <- dataset[ ,1:end_idx]
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)

splitted <- ratio_splitter(data, return.indicator = T)
trainset <- splitted$trainset
testset <- splitted$testset
train_indicator <- splitted$train_indicator
test_indicator <- splitted$test_indicator

result <- NULL 
num_factors <- 23 # number of latent factor for ageing data
param_grid <- expand.grid(lambda = seq(1, 12, length.out = 12), 
			  alpha = seq(0.1, 0.4, length.out = 4))

# for(i in seq(nrow(param_grid))){
#     cat('parameter grid:', paste(param_grid[i,], collapse = ','), "\n")
#     lambda <- round(param_grid[i, 1], 2)
#     alpha <- round(param_grid[i, 2], 2)
#     obj <- train_amf_model(trainset, confounders, num_factors = num_factors)
#     if(is.null(result)){
#         result <- c(round(param_grid[i, ], 2), obj$optimal_rmse)
#         result <- t(as.matrix(result))
#     }else{
#         result <- rbind(result, c(round(param_grid[i,], 2), obj$optimal_rmse))
#     }
#     write.csv(result, file = 'iMF_tuning_ageing_result.csv')
# }

# for(num_factors in seq(10, 30, by = 1)){
#     cat('length of factor vector:', num_factors, "\n")
#     obj <- train_amf_model(trainset, confounders, num_factors = num_factors)
# 
#     if(is.null(result)){
#         result <- c(num_factors, obj$optimal_rmse)
#         result <- t(as.matrix(result))
#     }else{
#         result <- rbind(result, c(num_factors, obj$optimal_rmse))
#     }
#     write.csv(result, file = 'iMF_tuning_num_factor_ageing.csv')
# }

optimal_rmse <- train_amf_model(trainset, confounders, num_factors = num_factors)
