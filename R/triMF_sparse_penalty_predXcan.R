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
source('./optimization_functions.R')
sourceCpp("./cdplusplus.cpp")
# source("cdplusplus.R")

# dump objects for easy debugging
dump_and_quit <- function() {
  # Save debugging info to file last.dump.rda
  dump.frames(to.file = TRUE)
  # Quit R with error status
  q(status = 1)
}
options(error = dump_and_quit)

predict <- function(disease_factor, tissue_factor, column_factor, cfd1, cfd2) {

    predictions <- (disease_factor[cfd1, ] + tissue_factor[cfd2, ] + disease_factor[cfd1, ] * tissue_factor[cfd2, ]) %*% column_factor
 
    return(predictions)
}

compute_loss <- function(disease_factor, tissue_factor, column_factor, sum_residual, verbose = T){

    # l2 penalty 
    row_reg_loss <- lambda * (norm(tissue_factor, "F")^2 + norm(disease_factor, "F")^2)
    col_reg_loss <- lambda * (1 - alpha) * norm(column_factor, "F")^2
    
    # l1 penalty
    l1_reg_loss <- lambda * alpha * sum(abs(column_factor))

    if(verbose){
        cat(paste0('total_residual:\t', sum_residual/2, 
                   '; \nrow_reg_loss:\t', row_reg_loss/2, 
                   '; \ncol_reg_loss:\t', col_reg_loss/2, 
                   '; \nl1_reg_loss:\t', l1_reg_loss, "\n"))
    }
    loss <- sum_residual/2 + row_reg_loss/2 + col_reg_loss/2 + l1_reg_loss
    return(loss)
}

evaluate_model <- function(trainset, testset, prediction, global_mean, iter = NA, verbose = T, debug = F){

    train_diff <- (trainset - prediction)
    sum_residual <- sum((train_diff[train_indicator])^2)
    train_rmse <- sqrt(sum_residual/sum(train_indicator))

    test_diff <- ((testset + global_mean - prediction)[test_indicator])^2
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

train_amf_model <- function(trainset, confounders, global_mean, num_factors = 30, num_iterations = 50000){
    # for disease tissue specific, confounders should have two columns, one for disease and another for tissue
    loss <- 0.0; last_loss <- Inf; rmse <- 0.0; optimal_rmse = Inf
    num_col <- ncol(trainset)
    col_loop <- seq(num_col)

    disease_factor <- matrix(generate_rand_vector(num_factors * num_disease), ncol = num_factors)
    tissue_factor <- matrix(generate_rand_vector(num_factors * num_tissue), ncol = num_factors)
    column_factor <- matrix(generate_rand_vector(num_factors * ncol(trainset)), nrow = num_factors)

    cfd1 <- as.integer(confounders[[1]])	
    cfd2 <- as.integer(confounders[[2]])

    predictions <- predict(disease_factor, tissue_factor, column_factor, cfd1, cfd2)
    metrics <- evaluate_model(trainset, testset, predictions, global_mean)
    loss <- compute_loss(disease_factor, tissue_factor, column_factor, metrics$sum_residual)

    for(iter in seq(num_iterations)){
        last_loss <- loss
        if(iter %% 10 == 0){
            cat("Iteration ", iter, '------------------------------------------\n')	
            cat('sum of squared error in iteration ', iter, ' before optimizating stage_factor: ', last_loss, '\n')
        }
        
        # fix all other parameters and update disease factor
        tissue_matrix <- tissue_factor %*% column_factor
        residual <- trainset - tissue_matrix[cfd2, ]
        residual[zero_idx] <- 0
        disease_factor <- optimize_row_l2(residual, disease_factor, tissue_factor, column_factor, cfd1, cfd2, lambda, n_cores = 16)
	    
        # fix all other parameters and update tissue factors
        disease_matrix <- disease_factor %*% column_factor
        residual <- trainset - disease_matrix[cfd1, ]
        residual[zero_idx] <- 0
        tissue_factor <- optimize_row_l2(residual, tissue_factor, disease_factor, column_factor, cfd2, cfd1, lambda, n_cores = 16)

        # fix all other parameters and update column factor term
        row_factor <- disease_factor[cfd1, ] + tissue_factor[cfd2, ] + disease_factor[cfd1, ] * tissue_factor[cfd2, ]
        column_factor <- optimize_col(trainset, row_factor, column_factor, lambda, alpha, n_cores = 30, tol = 1e-5 * decay)

        if(iter %% 10 == 0){
            predictions <- predict(disease_factor, tissue_factor, column_factor, cfd1, cfd2)
            metrics <- evaluate_model(trainset, testset, predictions, global_mean)
            loss <- compute_loss(disease_factor, tissue_factor, column_factor, metrics$sum_residual)

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
                        trainset = trainset,
                        disease_factor = disease_factor,
                        tissue_factor = tissue_factor,
                        column_factor = column_factor, 
                        optimal_rmse = optimal_rmse)

    return(fitted_obj)
}

load('~/data/zscore_matrix_selected.RData')

zscore[is.na(zscore)] <- 0
splitted <- ratio_splitter(zscore, return.indicator = T)
trainset <- splitted$trainset
testset <- splitted$testset
train_indicator <- (trainset != 0)
test_indicator <- (testset != 0)

names <- rownames(trainset)
grouping <- t(unname(sapply(names, function(s) split_str(s))))
# grouping <- grouping[grouping[ , 1] == single_disease, ]
confounders <- as.data.frame(grouping, stringsAsFactors = T)
colnames(confounders) <- c('disease', 'tissue')
num_disease <- nlevels(confounders[[1]])
num_tissue <- nlevels(confounders[[2]])
levels(confounders[[1]]) <- seq(num_disease)
levels(confounders[[2]]) <- seq(num_tissue)
confounders$disease <- as.integer(confounders[[1]])
confounders$tissue <- as.integer(confounders[[2]])

num_factors <- 50

global_mean <- calculate_grand_mean(trainset)
cat("Global mean:", global_mean, "\n")

# for(i in seq(nrow(param_grid))){
#     cat('parameter grid:', paste(param_grid[i,], collapse = ','), "\n")
#     lambda <- round(param_grid[i, 1], 2)
#     alpha <- round(param_grid[i, 2], 2)
#     obj <- train_amf_model(trainset - global_mean * train_indicator, confounders, global_mean, num_factors = num_factors)
#     if(is.null(result)){
#         result <- c(round(param_grid[i, ], 2), obj$optimal_rmse)
#         result <- t(as.matrix(result))
#     }else{
#         result <- rbind(result, c(round(param_grid[i,], 2), obj$optimal_rmse))
#     }
#     write.csv(result, file = 'iMF_tuning_predXcan_result.csv')
# }

for(num_factors in seq(20, 50, by = 1)){
    cat('length of factor vector:', num_factors, "\n")
    obj <- train_amf_model(trainset - global_mean * train_indicator, confounders, global_mean, num_factors = num_factors)

    if(is.null(result)){
        result <- c(num_factors, obj$optimal_rmse)
        result <- t(as.matrix(result))
    }else{
        result <- rbind(result, c(num_factors, obj$optimal_rmse))
    }
    write.csv(result, file = 'iMF_tuning_num_factor_predXcan.csv')
}

# optimal_rmse <- train_amf_model(trainset, confounders, num_factors = num_factors)

# fitted_obj <- train_amf_model(trainset - global_mean * train_indicator, confounders, global_mean, num_factors = num_factors)
