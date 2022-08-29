#' create an insider object with the provided parameters
#'
#' @param data the whole data for insider, including trainset and testset
#' @param confounder a matrix with each column indicating the belonging for each covariate
#' @param split_ratio a proportion of elements from data considered as testset
#' @param global_tol the global convergence criteria
#' @param sub_tol  the convergence creteria for each elastic net problem
#' @param tuning_iter number of steps will run in tuning
#' @param max_iter maxiumme number of steps will run in fitting insider
#'
#' @return an insider object
#' @export
#'
#' @examples ..
insider <- function(data, confounder, interaction_idx, split_ratio = 0.1, global_tol = 1e-9, sub_tol = 1e-5, tuning_iter = 30, max_iter = 50000){

    # split data into two pieces
    dataset <- ratio_splitter(data, ratio = split_ratio)
    trainset <- dataset[['trainset']]
    testset <- dataset[['testset']]

    # create insider class
    object <- structure(list(), class = "insider")
    object[['data']] <- data

    if(is.integer(interaction_idx) & (length(interaction_idx) > 1)){
        
        if(max(interaction_idx) > ncol(confounder)){
            stop("The interaction_idx is out of the range of confounder!")
        }

        unique_cfd <- unique(confounder[, interaction_idx])

        interaction_indicator <- rep(0, nrow(confounder))
        for(k in 1:nrow(unique_cfd)){
            selected <- apply(confounder[, interaction_idx], 1, function(x) all(x == unique_cfd[k,]))
            interaction_indicator[selected] <- k
        }
    }else{
        stop("The interaction_idx should be integers and its length must be greater than or equal to 2!")
    }
    
    object[['confounder']] <- cbind(confounder, interaction_indicator)

    # generate indicator for easy operation in C++
    object[['train_indicator']] <- apply(dataset[['train_indicator']], 2, as.integer)

    params <- list(global_tol = global_tol, sub_tol = sub_tol,
                   tuning_iter = tuning_iter, max_iter = max_iter)

    object[['params']] <- params

    return(object)
}

#' tune hyperparameters for insider object
#'
#' @param object an insider object
#' @param latent_dimension the rank of latent representations
#' @param lambda l2-regularization for insider
#' @param alpha l1-regularization for insider
#'
#' @return tuning results
#' @export
#'
#' @examples ..
tune <- function(object, latent_dimension = NULL, lambda = 1.0, alpha = 0.1){
    
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

            cat('Latent rank: ', latent_rank, "---------------------------------\n")

            # create a matrix for each covariate and put them into a List
            confounder_num <- ncol(object[['confounder']])
            confounder_list <- lapply(1:confounder_num, function(i){
                factor_num <- length(unique(object[['confounder']][,i]))
                matrix(init_parameters(factor_num * latent_rank), ncol = latent_rank)
            })
            
            column_factor <- matrix(init_parameters(latent_rank * ncol(data)), nrow = latent_rank)

            fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                    latent_rank, 1.0, 0.1, 1, global_tol, sub_tol, tuning_iter);

            if(is.null(rank_tuning)){
                rank_tuning <- c(latent_rank, fitted_obj$train_rmse, fitted_obj$test_rmse)
                rank_tuning <- t(as.matrix(rank_tuning))
            }else{
                rank_tuning <- rbind(rank_tuning, c(latent_rank, fitted_obj$train_rmse, fitted_obj$test_rmse))
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

        confounder_num <- ncol(object[['confounder']])
        confounder_list <- lapply(1:confounder_num, function(i){
            factor_num <- length(unique(object[['confounder']][,i]))
            matrix(init_parameters(factor_num * latent_rank), ncol = latent_rank)
        })
    
        column_factor <- matrix(init_parameters(latent_rank * ncol(data)), nrow = latent_rank)

        param_grid <- expand.grid(lambda = lambda, alpha = alpha)

        for(i in seq(nrow(param_grid))){
            cat('parameter grid:', paste(round(param_grid[i,], 2), collapse = ','), "---------------------------------\n")
            
            lambda <- round(param_grid[i, 1], 2)
            alpha <- round(param_grid[i, 2], 2)

            fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                                   latent_rank, lambda, alpha, 1, global_tol, sub_tol, tuning_iter);

            if(is.null(reg_tuning)){
                reg_tuning <- c(round(param_grid[i, ], 2), fitted_obj$train_rmse, fitted_obj$test_rmse)
                reg_tuning <- t(as.matrix(reg_tuning))
            }else{
                reg_tuning <- rbind(reg_tuning, c(round(param_grid[i,], 2), fitted_obj$train_rmse, fitted_obj$test_rmse))
            }
            write.csv(reg_tuning, file = 'insider_R', latent_rank, '_reg_tuning_result.csv')
        }
    }
    return(list(rank_tuning = rank_tuning, latent_rank = latent_rank, reg_tuning = reg_tuning))
}

#' fit insider object
#'
#' @param object an insider object
#' @param latent_dimension the rank of latent representations
#' @param lambda l2-regularization for insider
#' @param alpha l1-regularization for insider
#'
#' @return a fitted insider object
#' @export
#'
#' @examples ..
fit <- function(object, latent_dimension = NULL, lambda = NULL, alpha = NULL){
    
    global_tol <- object[['params']][['global_tol']]
    sub_tol <- object[['params']][['sub_tol']]
    max_iter <- object[['params']][['max_iter']]
    
    confounder_num <- ncol(object[['confounder']])
    confounder_list <- lapply(1:confounder_num, function(i){
        factor_num <- length(unique(object[['confounder']][,i]))
        matrix(init_parameters(factor_num * latent_dimension), ncol = latent_dimension)
    })

    column_factor <- matrix(init_parameters(latent_dimension * ncol(data)), nrow = latent_dimension)

    fitted_obj <- optimize(object[['data']], confounder_list, column_factor, object[['confounder']], object[['train_indicator']], 
                           latent_dimension, lambda, alpha, 0, global_tol, sub_tol, max_iter);

    object[['cfd_matrices']] <- fitted_obj[['row_matrices']]
    object[['column_factor']] <- fitted_obj[['column_factor']]
    object[['test_rmse']] <- fitted_obj[['test_rmse']]

    return(object)
}

#' @export 
capture_interaction <- function(object, inc_cfd, latent_dimension, lambda, alpha, ratio = 0.2){
    
    tuning <- ifelse((length(lambda) > 1) | (length(alpha) > 1), 1, 0)

    if(tuning == 1){
        dataset <- ratio_splitter(object[['data']], ratio = ratio)
        object[['train_indicator']] <- apply(dataset[['train_indicator']], 2, as.integer)
    }
    
    confounder_num <- ncol(object[['confounder']])
    residual <- object[['data']]

    for(i in 1:confounder_num){
        sub_predictions <- object[['cfd_matrices']][[i]] %*% object[['column_factor']]
        residual <- residual - sub_predictions[object[['confounder']][,i], ]
    }

    confounder <- object[['confounder']][, inc_cfd]
    unique_cfd <- unique(confounder)

    interaction_indicator <- rep(0, nrow(confounder))
    for(k in 1:nrow(unique_cfd)){
        selected <- apply(confounder, 1, function(x) all(x == unique_cfd[k,]))
        interaction_indicator[selected] <- k
    }


    rmse <- NULL
    param_grid <- expand.grid(lambda = lambda, alpha = alpha)

    if(tuning == 1){
        for(i in seq(nrow(param_grid))){

            cat('parameter grid:', paste(round(param_grid[i,], 2), collapse = ','), "---------------------------------\n")
                
            lambda <- round(param_grid[i, 1], 2)
            alpha <- round(param_grid[i, 2], 2)

            row_interaction <- matrix(0, nrow = nrow(confounder), ncol = latent_dimension)

            interactions <- fit_interaction(residual, object[['train_indicator']], interaction_indicator, object[['column_factor']], lambda, alpha, tuning)

            predictions <- interactions %*% object[['column_factor']]
            diff <- residual - predictions[interaction_indicator,]

            train_rmse <- sqrt(mean(diff[object[['train_indicator']] == 1]^2))
            test_rmse <- sqrt(mean(diff[object[['train_indicator']] == 0]^2))

            cat('train_rmse:', train_rmse, ", test_rmse:", test_rmse, '\n')

            if(is.null(rmse)){
                rmse <- c(lambda, alpha, train_rmse, test_rmse)
                rmse <- t(as.matrix(rmse))
            }else{
                rmse <- rbind(rmse, c(lambda, alpha, train_rmse, test_rmse))
            }
            write.csv(rmse, file = 'insider_interaction_tuning_result.csv')
        }
        rmse <- as.data.frame(rmse)
        colnames(rmse) <- c('lambda', 'alpha', 'train_rmse', 'test_rmse')
        
        min_idx <- which.min(rmse[['test_rmse']])
        lambda <- rmse[['lambda']][min_idx]
        alpha <- rmse[['alpha']][min_idx]
        tuning <- 0
    }

    cat('lambda: ', lambda, '; alpha: ', alpha, '\n')

    interactions <- fit_interaction(residual, object[['train_indicator']], interaction_indicator, object[['column_factor']], lambda, alpha, tuning)

    object[['interactions']] <- interactions
    return(object)
}


