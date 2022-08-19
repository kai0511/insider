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
# source('./global_conf.R')
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

predict <- function(stage_factor, tissue_factor, column_factor, cfd1, cfd2){
    predictions <- (stage_factor[cfd1, ] + tissue_factor[cfd2, ] + stage_factor[cfd1, ] * tissue_factor[cfd2, ]) %*% column_factor
    return(predictions)
}

compute_loss <- function(stage_factor, tissue_factor, column_factor, sum_residual, verbose = T){
    # l1 penalty
    l1_reg_loss <- alpha * lambda * sum(abs(column_factor))
   
    # l2 penalty 
    col_reg_loss <- (1 - alpha) * lambda * norm(column_factor, "F")^2
    row_reg_loss <- lambda * (norm(stage_factor, "F")^2 + norm(tissue_factor, "F")^2)

    if(verbose){
        cat(paste0("total_residual:\t", sum_residual/2, 
                   "; \nrow_reg_loss:\t", row_reg_loss/2, 
                   "; \ncol_reg_loss:\t", col_reg_loss/2, 
                   "; \nl1_reg_loss:\t", l1_reg_loss, "\n"))
    }
    loss <- sum_residual/2 + row_reg_loss/2 + col_reg_loss/2 + l1_reg_loss
    return(loss)
}

evaluate_model <- function(trainset, testset, prediction, iter = NA, verbose = T){

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

train_mf_model <- function(trainset, testset, confounders, num_factors = 50, num_iterations = 50000, verbose = T, debug = T){
    # for stage tissue specific, confounders should have two columns, one for stage and another for tissue
    loss <- 0.0; last_loss <- Inf; rmse <- 0.0; optimal_rmse = Inf; decay = 1
    # gradients <- 0.0; last_gradients <- Inf;

    num_col <- ncol(trainset)
    col_loop <- seq(num_col)
    zero_idx <- (trainset == 0)

    cfd1 <- as.integer(confounders[[1]])	
    cfd2 <- as.integer(confounders[[2]])
    num_stage <- length(unique(cfd1))
    num_tissue <- length(unique(cfd2))

    stage_factor <- matrix(generate_rand_vector(num_factors * num_stage), ncol = num_factors)
    tissue_factor <- matrix(generate_rand_vector(num_factors * num_tissue), ncol = num_factors)
    column_factor <- matrix(generate_rand_vector(num_factors * ncol(trainset)), nrow = num_factors)

    indicator <- apply(train_indicator, 2, as.numeric)
    predictions <- predict(stage_factor, tissue_factor, column_factor, cfd1, cfd2)
    metrics <- evaluate_model(trainset, testset, predictions)
    loss <- compute_loss(stage_factor, tissue_factor, column_factor, metrics$sum_residual)

    for(iter in seq(num_iterations)){
        last_loss <- loss
        if(iter %% 10 == 0){
            cat("Iteration ", iter, '------------------------------------------\n')	
            cat('sum of squared error in iteration ', iter, ' before optimizating stage_factor: ', last_loss, '\n')
        }
        
        # fix all other parameters and update stage factor
        # tissue_matrix <- tcrossprod(tissue_factor, column_factor)
        # stage_factor_update <- mclapply(unique(cfd1), function(d) {

        #     idx <- which(cfd1 == d)
        #     ord <- cfd2[idx]
        #     mat <- trainset[idx, , drop = F] - tissue_matrix[ord, , drop = F]
        #     num_row <- nrow(mat)
        #     indication_mat <- train_indicator[idx, , drop = F]

        #     M_matrix <- lapply(seq(num_row), function(i){
        #         sweep(column_factor[indication_mat[i, ],], MARGIN=2, 1 + tissue_factor[ord[i], ], '*')
        #     })
        #     left_matrix <- as.matrix(do.call(rbind, M_matrix))
        #     right_vector <- as.vector(t(mat))[as.vector(t(indication_mat))]

        #     d_factor <- chol2inv(chol(crossprod(left_matrix) + diag(r1_lambda, num_factors, num_factors))) %*% crossprod(left_matrix, right_vector)
        #     as.vector(d_factor)
        # }, mc.cores = 10)
        # stage_factor <- do.call(rbind, stage_factor_update)
        tissue_matrix <- tissue_factor %*% column_factor
        residual <- trainset - tissue_matrix[cfd2, ]
        # residual[zero_idx] <- 0
        stage_factor <- optimize_row_l2(residual, indicator, stage_factor, tissue_factor, column_factor, cfd1, cfd2, lambda, n_cores = 16)
        # stage_factor <- result[["updates"]]
        # stage_gradient <- result[["gradients"]]

        # fix all other parameters and update tissue factors
        # stage_matrix <- tcrossprod(stage_factor, column_factor)
        # tissue_factor_update <- mclapply(unique(cfd2),  function(t){

        #     idx <- which(cfd2 == t)
        #     ord <- cfd1[idx]
        #     mat <- trainset[idx, , drop = F] - stage_matrix[ord, , drop = F]
        #     num_row <- nrow(mat)
        #     indication_mat <- train_indicator[idx, , drop = F]

        #     M_matrix <- lapply(seq(num_row), function(i){
        #         sweep(column_factor[indication_mat[i, ],], MARGIN=2, 1 + stage_factor[ord[i], ], '*')
        #     })
        #     left_matrix <- do.call(rbind, M_matrix)
        #     right_vector <- as.vector(t(mat))[as.vector(t(indication_mat))]

        #     # coordiante descent part
        #     t_factor <- chol2inv(chol(crossprod(left_matrix) + diag(r1_lambda, num_factors, num_factors))) %*% crossprod(left_matrix, right_vector)
        #     as.vector(t_factor)
        # }, mc.cores = 10)
        # tissue_factor <- do.call(rbind, tissue_factor_update)
        stage_matrix <- stage_factor %*% column_factor
        residual <- trainset - stage_matrix[cfd1, ]
        # residual[zero_idx] <- 0

        tissue_factor <- optimize_row_l2(residual, indicator, tissue_factor, stage_factor, column_factor, cfd2, cfd1, lambda, n_cores = 16)
        # tissue_factor <- result[["updates"]]
        # tissue_gradient <- result[["gradients"]]

        # fix all other parameters and update column factors
        row_factor <- stage_factor[cfd1, ] + tissue_factor[cfd2, ] + stage_factor[cfd1, ] * tissue_factor[cfd2, ]
        column_factor <- optimize_col(trainset, indicator, row_factor, column_factor, lambda, alpha, n_cores = 30, tol = 1e-5 * decay)
        # column_factor <- result[["updates"]]
        # column_gradient <- result[["gradients"]]
        # column_gradient[column_factor == 0] <- 0

        # gradients <- sum(abs(column_gradient)) + sum(abs(tissue_gradient)) + sum(abs(stage_gradient))
        # cat(paste0('sum of absolute gradients in iteration ', iter, ' after optimization: ', gradients, '\n'))
        # cat(paste0('delta gradients in iteration ', iter, ' after optimization: ', last_gradients - gradients, '\n'))
        # last_gradients <- gradients 

        # calculate residual

        if(iter %% 10 == 0){
            last_loss <- loss
	    predictions <- predict(stage_factor, tissue_factor, column_factor, cfd1, cfd2)
            metrics <- evaluate_model(trainset, testset, predictions, iter = iter)
            loss <- compute_loss(stage_factor, tissue_factor, column_factor, metrics$sum_residual)

            rmse <- metrics[['test_rmse']]
            cat(paste0('sum of squared error in iteration ', iter, ' after optimizing factor: ', rmse, '\n'))
            
            decay <- ifelse((last_loss - loss)/100 >= 1, 1, signif((last_loss - loss)/100, digits=1))
            if(is_converged(loss, last_loss, iter, learner = 'ibiasedMF_L1_penalty', thres = 1e-10)) {
                break
            }
            optimal_rmse <- rmse
        }
        
    }

    fitted_obj <- list (iter = iter,
                        predictions = predictions,
                        stage_factor = stage_factor,
                        tissue_factor = tissue_factor,
                        column_factor = column_factor, 
                        optimal_rmse = optimal_rmse)

    save(fitted_obj, file = paste0('gene_expression_iMF_L1L2_penalty_',  num_factors, '_.RData'))
    return(fitted_obj)
}



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

splitted <- ratio_splitter(data)
trainset <- splitted$trainset
testset <- splitted$testset
train_indicator <- splitted$train_indicator
test_indicator <- splitted$test_indicator

result <- NULL 
num_factors <- 32  # number of latent factor for brainspan data
# parameter grid for brainspan dataset
param_grid <- expand.grid(lambda = seq(10, 30, length.out = 10), 
			  # alpha = seq(0, 1, length.out = 11))
			  alpha = c(0.5, 0.6, 0.7, 0.8))

for(i in seq(nrow(param_grid))){
    cat('parameter grid:', paste(param_grid[i,], collapse = ','), "\n")
    lambda <- round(param_grid[i, 1], 2)
    alpha <- round(param_grid[i, 2], 2)
    obj <- train_mf_model(trainset, testset, confounders, num_factors = num_factors)
    if(is.null(result)){
        result <- c(round(param_grid[i, ], 2), obj$optimal_rmse)
        result <- t(as.matrix(result))
    }else{
        result <- rbind(result, c(round(param_grid[i,], 2), obj$optimal_rmse))
    }
    write.csv(result, file = 'iMF_tuning_brainspan_result.csv')
}

# for(num_factors in seq(15, 40, by = 1)){
#     cat('length of factor vector:', num_factors, "\n")
#     obj <- train_mf_model(trainset, testset, confounders, num_factors = num_factors)
# 
#     if(is.null(result)){
#         result <- c(num_factors, obj$optimal_rmse)
#         result <- t(as.matrix(result))
#     }else{
#         result <- rbind(result, c(num_factors, obj$optimal_rmse))
#     }
#     write.csv(result, file = 'iMF_tuning_num_factor_brainspan.csv')
# }

# obj <- train_mf_model(trainset, testset, confounders, num_factors = num_factors)
