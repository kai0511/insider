insider <- function(data, confounder, split_ratio = 0.1, global_tol = 1e-9, sub_tol = 1e-5, tuning_iter = 100, max_iter = 50000){

    # split data into two pieces
    dataset <- ratio_splitter(data, ratio = split_ratio)
    trainset <- dataset[['trainset']]
    testset <- dataset[['testset']]

    # generate indicator for easy operation in C++
    train_indicator <- apply(dataset[['train_indicator']], 2, as.numeric)

    # create insider class
    object <- structure(list(), class = "insider")
    object[['data']] <- data
    object[['confounder']] <- confounder
    object[['trainset']] <- trainset
    object[['testset']] <- testset
    object[['train_indicator']] <- train_indicator

    params <- list(global_tol = global_tol, sub_tol = sub_tol,
                   tuning_iter = tuning_iter, max_iter = max_iter)

    object[['params']] <- params

    return(object)
}

tune.insider <- function(object, latent_dimension = NULL, lambda = 1.0, alpha = 0.1){
    
    if(!is.integer(latent_dimension) | !is.numeric(lambda) | !is.numeric(alpha)){
        stop("TUNNING: The element of latent_dimension, lambda, and alpha should be integer, numeric, and numeric.")
    }

    if((length(latent_dimension) <= 1) & (length(lambda) <= 1 | length(alpha) <= 1)){
        stop("TUNNING: The length of either latent_dimension or lambda and alpha should be greater than 1.")
    }

    global_tol <- object[['params']][['global_tol']]
    sub_tol <- object[['params']][['sub_tol']]
    tuning_iter <- object[['params']][['tuning_iter']]

    rank_tuning <- NULL; reg_tuning <- NULL

    # tune the rank of latent representations
    if(length(latent_dimension) > 1){

        for(latent_rank in latent_dimension){
            # create a matrix for each covariate and put them into a List
            confounder_num <- ncol(confounder)
            confounder_list <- lapply(1:confounder_num, function(i){
                factor_num <- unique(confounder[,i])
                matrix(init_parameters(factor_num * latent_rank), ncol = latent_rank)
            })
            
            column_factor <- matrix(init_parameters(latent_rank * ncol(data)), nrow = latent_rank)

            fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                    latent_rank, 1.0, 0.1, 1, global_tol, sub_tol, tuning_iter);

            if(is.null(rank_tuning)){
                rank_tuning <- c(latent_rank, fitted_obj$optimal_rmse)
                rank_tuning <- t(as.matrix(rank_tuning))
            }else{
                rank_tuning <- rbind(rank_tuning, c(latent_rank, fitted_obj$optimal_rmse))
            }
            write.csv(rank_tuning, file = 'insider_rank_tuning_result.csv')
        }
    }

    # select the latent rank with the lowest optimal rmse
    if(length(latent_dimension) > 1){
        latent_rank <- latent_dimension[which.min(rank_tuning[,2])]
    }else{
        latent_rank <- latent_dimension
    }
    
    # tune lambda and alpha with the selected latent rank
    if(length(lambda) > 1 | length(alpha) > 1){

        confounder_num <- ncol(confounder)
        confounder_list <- lapply(1:confounder_num, function(i){
            factor_num <- unique(confounder[,i])
            matrix(init_parameters(factor_num * latent_rank), ncol = latent_rank)
        })
    
        column_factor <- matrix(init_parameters(latent_rank * ncol(data)), nrow = latent_rank)

        param_grid <- expand.grid(lambda = lambda, alpha = alpha)

        for(i in seq(nrow(param_grid))){
            cat('parameter grid:', paste(round(param_grid[i,]), collapse = ','), "\n")
            
            lambda <- round(param_grid[i, 1], 2)
            alpha <- round(param_grid[i, 2], 2)

            fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                                   latent_rank, lambda, alpha, 1, global_tol, sub_tol, tuning_iter);

            if(is.null(reg_tuning)){
                reg_tuning <- c(round(param_grid[i, ], 2), fitted_obj$optimal_rmse)
                reg_tuning <- t(as.matrix(reg_tuning))
            }else{
                reg_tuning <- rbind(reg_tuning, c(round(param_grid[i,], 2), fitted_obj$optimal_rmse))
            }
            write.csv(reg_tuning, file = 'insder_reg_tuning_result.csv')
        }
    }
    return(list(rank_tuning = rank_tuning, latent_rank = latent_rank, reg_tuning = reg_tuning))
}

fit.insider <- function(object, latent_dimension = NULL, lambda = NULL, alpha = NULL){
    
    global_tol <- object[['params']][['global_tol']]
    sub_tol <- object[['params']][['sub_tol']]
    max_iter <- object[['params']][['max_iter']]

    confounder_num <- ncol(confounder)
    confounder_list <- lapply(1:confounder_num, function(i){
        factor_num <- unique(confounder[,i])
        matrix(init_parameters(factor_num * latent_rank), ncol = latent_rank)
    })

    column_factor <- matrix(init_parameters(latent_rank * ncol(data)), nrow = latent_rank)

    fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                           latent_rank, lambda, alpha, 1, global_tol, sub_tol, max_iter);

    object[['cfd_matrices']] <- fitted_obj[['row_matrices']]
    object[['column_factor']] <- fitted_obj[['column_factor']]
    object[['optimal_rmse']] <- fitted_obj[['optimal_rmse']]

    return(object)
}