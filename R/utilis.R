require(methods)

`%+%` <- function(m, v){
    # add vector v to matrix m by column
    if(class(v) == 'matrix'){
        v <- as.vector(v)
    }  
    sweep(m, 2, v, '+')
}

`%-%` <- function(m, v){
    # minus vector v from matrix m by column
    if(class(v) == 'matrix'){
        v <- as.vector(v)
    } 
    sweep(m, 2, v, '-')
}

calculate_grand_mean <- function(data){
    # here data should be a R matrix 
    # return(mean(data, na.rm = T))
    return(mean(data[data != 0]))
}

truncated_sim <- function(x){
    lower_truncated <- pmax(x, -1)
    truncated_x <- pmin(lower_truncated, 1)
    return(truncated_x)
}

find_top_n_indices <- function(x, n){
    # this method is Deprecated
    # set the quantile to correspond to n top elements
    quant <- n / length(x)
    #select the cutpoint to get the quantile above quant
    cp <- quantile(x, probs = 1.0 - quant, na.rm = T)
    #select the elements above the cutpoint
    return(which(x > cp[1]))
}

calculate_idx <- function(idx, num_row){
    # calculate the position in trainset matrix
    col_idx <- floor(idx / num_row)  # column idx
                
    if(idx == col_idx * num_row){
        row_idx <- num_row
    } else {
        row_idx <- idx - col_idx * num_row
        col_idx <- col_idx + 1
    }
    return(c(row_idx, col_idx))
}

generate_rand_vector <- function(size, init_mean = 0.0, init_std = 0.001) {
    # init_mean = 0.0; init_std <- 0.001
    return(rnorm(size, mean = init_mean, sd = init_std))
}

split_str <- function(s){
    # l <- unlist(strsplit(s, split=c('-|_')))
    l <- unlist(strsplit(s, split=c('_')))
    idx <- which(l == 'v7')
    len <- length(l)
    
    disease <- l[1]
    tissue <- paste(l[(idx + 1): len], collapse = '_')
    c(disease, tissue)
}

obtain_indication_matrix <- function(trainset, only_positive = F){
    indication_matrix <- matrix(0, nrow = nrow(trainset), ncol = ncol(trainset))
    if(only_positive){
        indication_matrix[!is.na(trainset)] <- 1
    } else {
        indication_matrix[!is.na(trainset)] <- 1
        indication_matrix[trainset < 0] <- -1
    }
    return(indication_matrix)
}

construct_sym_matrix <- function(dimension, row_names = NULL, col_names = NULL){
    tmp_matrix <- matrix(generate_bias_vector(dimension * dimension), ncol = dimension)
    sym_matrix <- (tmp_matrix + t(tmp_matrix))/2
    colnames(sym_matrix) <- row_names
    rownames(sym_matrix) <- col_names
    return(sym_matrix)
}

ratio_splitter <-function(data, ratio = 0.1, rm.na.col = T){
    # default ratio for test is 0.1, that is, 10% obs. will be randomly assigned to testset
    
    train_indicator <- matrix(T, nrow = nrow(data), ncol = ncol(data))
    testset <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    test_indicator <- matrix(F, nrow = nrow(data), ncol = ncol(data))

    test_idx <- sample(seq(length(data)), floor(length(data) * ratio), replace = F)
    testset[test_idx] <- data[test_idx]
    data[test_idx] <- 0

    train_indicator[test_idx] <- F
    # train_indicator <- as(train_indicator, "sparseMatrix")
    test_indicator[test_idx] <- T
    # test_indicator <- as(test_indicator, "sparseMatrix")

    num_per_col <- apply(data, 2, function(x) sum(x != 0))
    cat(paste0('number of all zero columns removed: ',  sum(num_per_col == 0)), '\n')
    if (rm.na.col) {
        return(list(trainset = data[,num_per_col != 0], 
                    testset = testset[, num_per_col != 0],
                    train_indicator = train_indicator[,num_per_col != 0], 
                    test_indicator = test_indicator[,num_per_col != 0]))
    } else{
        return(list(trainset = data, 
                    testset = testset,
                    train_indicator = train_indicator, 
                    test_indicator = test_indicator))
    }
}

update_learn_rate <- function(loss, last_loss, learn_rate, iter) {
    
    if (is_bold_driver & iter > 1) {
        learn_rate <- ifelse(abs(last_loss) > abs(loss), learn_rate * 1.05, learn_rate * 0.9)
    } else if (decay > 0 & decay < 1) {
        learn_rate <- learn_rate * decay
    }

    # limit to max_learn_rate after update
    if (max_learn_rate > 0 & learn_rate > max_learn_rate) {
        learn_rate <- max_learn_rate
    }
    return(learn_rate)
}

update_lr <- function(loss, last_loss, learn_rate, learn_rate1, iter) {
    
    if (is_bold_driver & iter > 1) {
        learn_rate <- ifelse(abs(last_loss) > abs(loss), learn_rate * 1.05, learn_rate * 0.9)
        learn_rate1 <- ifelse(abs(last_loss) > abs(loss), learn_rate1 * 1.05, learn_rate1 * 0.9)
    } else if (decay > 0 & decay < 1) {
        learn_rate <- learn_rate * decay
        learn_rate1 <- learn_rate1 * decay
    }

    # limit to max_learn_rate after update
    learn_rate <- ifelse(learn_rate > max_learn_rate, max_learn_rate, learn_rate)
    learn_rate1 <- ifelse(learn_rate1 > max_learn_rate, max_learn_rate, learn_rate1)
    return(c(learn_rate, learn_rate1))
}

is_converged <- function(loss, last_loss, iter, learner = 'biasedMF', thres = 1e-8, verbose = T){
    delta_loss <- last_loss - loss

    if (verbose) {
        info <- paste0(learner, " iter ", iter,  ": loss = ", loss, ", delta_loss = ", delta_loss)
        cat(info, '\n')
    }

    if(is.na(loss) | is.infinite(loss)){
        cat("Loss = NaN or Infinity: current settings does not fit! Change the settings and try again!")
    }
    return (abs(delta_loss)/loss < thres)
}

compute_interaction_loss <- function(stage_vec, tissue_vec, column_vec, sum_residual, verbose = T){
    l2_reg <- r2_lambda * (sum(stage_vec^2) + sum(tissue_vec^2) + sum(column_vec^2))

    if(verbose){
        cat(paste0("total_residual:\t", sum_residual/2, 
                   "; \nl2_reg:\t", l2_reg/2, "\n"))
    }

    loss <- sum_residual/2 + l2_reg/2
    return(loss)
}


fit_interaction <- function(trainset, testset, prediction, confounders, num_iterations = 20000, verbose = T){
    # for stage tissue specific, confounders should have two columns, one for stage and another for tissue
    loss <- 0.0; last_loss <- Inf; rmse <- 0.0; optimal_rmse = Inf
    num_col <- ncol(residual)
    col_loop <- seq(num_col)

    train_zero_idx <- (trainset == 0)
    train_residual <- trainset - prediction
    train_residual[train_zero_idx] <- 0

    test_zero_idx <- (testset == 0)
    test_residual <- testset - prediction
    test_residual[test_zero_idx] <- 0
    train_obs <- sum(train_indicator)

    cfd1 <- as.integer(confounders[[1]])	
    cfd2 <- as.integer(confounders[[2]])
    num_stage <- length(unique(cfd1))
    num_tissue <- length(unique(cfd2))

    stage_vec <- generate_rand_vector(num_stage)
    tissue_vec <- generate_rand_vector(num_tissue)
    column_vec <- generate_rand_vector(ncol(residual))

    interaction <- tcrossprod(stage_vec[cfd1] * tissue_vec[cfd2], column_vec)
    metrics <- evaluate_model(trainset, testset, interaction)
    loss <- compute_interaction_loss(stage_vec, tissue_vec, column_vec, metrics$sum_residual)

    for(iter in seq(num_iterations)){
        
        stage_vec <- optimize_row_interaction(train_residual, stage_vec, cfd1, tissue_vec, cfd2, column_vec, r2_lambda)
        tissue_vec <- optimize_row_interaction(train_residual, tissue_vec, cfd2, stage_vec, cfd1, column_vec, r2_lambda)
    
        row_interaction <- stage_vec[cfd1] * tissue_vec[cfd2]
        column_vec <- optimize_col_interaction(train_residual, row_interaction, column_vec, r2_lambda)

        interaction <- tcrossprod(stage_vec[cfd1] * tissue_vec[cfd2], column_vec)

        metrics <- evaluate_model(trainset, testset, interaction)
        rmse <- metrics[['test_rmse']]
        if(optimal_rmse > rmse) optimal_rmse <- rmse
        
        last_loss <- loss
        loss <- compute_interaction_loss(stage_vec, tissue_vec, column_vec, metrics$sum_residual)

        if(is_converged(loss, last_loss, iter, learner = 'interaction model', thres = 1e-5)) {
            break
        }
    }
    obj <- list(iter = iter, 
                interaction = interaction, 
                stage_vec = stage_vec, 
                tissue_vec = tissue_vec, 
                column_vec = column_vec, 
                optimal_rmse = optimal_rmse)
    
    return(obj)
}
