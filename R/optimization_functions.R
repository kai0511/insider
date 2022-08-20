
prox_l1 <- function(x, lambda){
    # proximal operator for lasso
    pmax(0, x - lambda) - pmax(0, -x - lambda)
}

proximal_gradient <- function(A, b, lambda, alpha, ABSTOL = 1e-8, MAX_ITER = 100){
    # reference see http://www.princeton.edu/~yc5/ele538b_sparsity/lectures/lasso_algorithm_extension.pdf
    # https://web.stanford.edu/~boyd/papers/prox_algs/lasso.html
    mu <-  1  
    beta <- 0.5

    AtA <- crossprod(A)
    Atb <- crossprod(A, b)
    x <- rep(0, ncol(A))
    x_prev <- x
    loss_prev <- Inf

    for(k in 1:MAX_ITER){
        y <- x + (k/(k+3))*(x - x_prev);
        repeat{
            grad_y <- drop(AtA %*% y - Atb)
            z <- 1/(1 + mu * lambda * (1 - alpha)) * prox_l1(y - mu * grad_y, mu * lambda * alpha)
            ssl <- 0.5 * sum((A %*% y - b)^2) + grad_y %*% (z - y) + (1/(2*mu))* sum((z - y)^2)
            if(0.5 * sum((A %*% z - b)^2) <= ssl){
                break
            }
            mu <- beta * mu;
        }
        x_prev <- x
        x <- z

        loss <- objective(A, b, x, lambda, alpha)
        if(k > 1 & abs(loss_prev - loss) < ABSTOL){
            break
        }
        loss_prev <- loss
    }
    return(x)
}

soft_threshold <- function(v){
    sapply(v, function(i) max(i, 0))
}

objective <- function(A, b, beta, lambda, alpha){
    # the loss for the following objective:
    # Fomular: 1/2 (A * beta - b)^2 + \lambda [(1-\alpha) ||beta||_2^2 + \alpha ||beta||_1]
    loss <- sum((A %*% beta - b)^2)/2 + lambda * ((1 - alpha) * sum(beta^2)/2 + alpha * sum(abs(beta)))
    return(loss)
}

compute_sub_loss <- function(residual, x, lambda, alpha, verbose = F){
    total_residual <- 0.5 * sum(residual^2)
    l2_loss <- 0.5 * (1 - alpha) * lambda * sum(x^2)
    l1_loss <- alpha * lambda * sum(abs(x))
    loss <- total_residual + l2_loss + l1_loss

    if(verbose){
        cat(paste0('total_residual:\t', total_residual, 
                   '; \nl2_loss:\t', l2_loss, 
                   '; \nl1_loss:\t', l1_loss, "\n"))
    }
    return(loss)
}

coordinate_descent <- function(X, y, lambda, alpha, beta, ex_idx, thres = 1e-10, verbose = FALSE){
    # for detail of optimization algorithm see: https://ieeexplore.ieee.org/document/6413853
    # screening rules to hasten computation. See: https://statweb.stanford.edu/~tibs/ftp/strong.pdf
    XtX <- crossprod(X)
    beta[ex_idx] <- 0
    residual <- drop(y - X %*% beta)
    pre_loss <- compute_sub_loss(residual, beta, lambda, alpha)

    while(1){
        for(k in sample.int(ncol(X))){
            if(k %in% ex_idx){
                next
            }
            
            upper_scala <- residual %*% X[, k] + beta[k] * XtX[k, k]
            right_scala <- XtX[k, k] + lambda * (1 - alpha)
            
            # soft-thresholding
            update <- sign(upper_scala) * max(abs(upper_scala) - lambda * alpha, 0) / right_scala
            residual <- residual - (drop(update) - beta[k]) * X[, k]
            beta[k] <- update
        }

        
        iter_loss <- compute_sub_loss(residual, beta, lambda, alpha)
        if(verbose) cat('delata loss: ', abs(iter_loss - pre_loss), '\n')

        # check convergence
        if(abs(iter_loss - pre_loss) < thres) break
        pre_loss <- iter_loss
    }
    return(beta)
}

safe_cd <- function(X, y, lambda, alpha, thres = 1e-10){
    # for detail of optimization algorithm see: https://ieeexplore.ieee.org/document/6413853
    # screening rules to hasten computation. See: https://statweb.stanford.edu/~tibs/ftp/strong.pdf
    proj <- abs(crossprod(X, y))
    XtX <- crossprod(X)
    max_lambda <- max(proj)
        
    excluded <- (proj < alpha * (2 * lambda - max_lambda))
    ex_idx <- seq(ncol(X))[excluded]
    
    # warm start
    beta <- drop(chol2inv(chol(XtX + (1 - alpha) * lambda * diag(1, ncol(X)))) %*% crossprod(X, y))
    beta[excluded] <- 0
    
    
    while(1){
        beta <- cd(X, y, lambda, alpha, beta, ex_idx, thres = thres)

        zero_indx <- which(beta == 0)
        lft <- crossprod(X, y - X %*% beta) - lambda * (1 - alpha) * beta
        rght <- sign(beta) * alpha * lambda
        rght[zero_indx] <- alpha * lambda

        # the key line to find the variable violating KKT condition
        idx <- (abs(lft[zero_indx]) > rght[zero_indx])

        if(sum(idx) != 0){
            ex_idx <- setdiff(ex_idx, zero_indx[idx])
        }else{
            break 
        }
    }
    return(beta)
}

feature_sign_with_screening <- function(X, y, sugg_start, lambda, alpha, 
                                        XtX, Xty, max_iter = 1000){
    #  update (20220610)
    #     adopt screening rules to filter out coefficients that shrink to zeros.
    #     for details, see https://statweb.stanford.edu/~tibs/ftp/strong.pdf. 
    #     use analytical solution to the elastic net.
    
    iter <- 0; col_num <- ncol(X)
    optimality1 <- FALSE; optimality0 <- FALSE
    pre_loss <- 0.0; inner_loss <- 0.0; line_search <- 0.0
    
    beta <- sugg_start; theta <- sign(beta)
    active_set <- rep(1, col_num)

    cutoff <- alpha * (2 * lambda - max(abs(Xty)))
    exc_idx <- (abs(Xty) < cutoff)
    active_set[exc_idx] <- 0
    beta[exc_idx] <- 0

    while(iter < max_iter){

        cat("Iter:", iter, "\n")
        
        while(TRUE) {
            line_search = 0.0
            inc_idx <- seq(col_num)[active_set == 1]

            X2 <- X[, inc_idx]
            inv_gram <- chol2inv(chol(XtX[inc_idx, inc_idx] + (1 - alpha) * lambda * diag(1, ncol(X2))))
            
            new_beta <- as.vector(inv_gram %*% (Xty[inc_idx] - lambda * alpha * theta[inc_idx]))

            cat("positions with different sign:", which(sign(new_beta) != theta[inc_idx]), "\n")
            idx <- (sign(new_beta) != theta[inc_idx])
            cat("pre_beta:", beta[inc_idx[idx]], "\n")
            cat("new_beta:", new_beta[sign(new_beta) != theta[inc_idx]], "\n")


            if (all(sign(new_beta) == theta[inc_idx])){
                beta[inc_idx] = new_beta
                optimality1 = TRUE
                line_search = 1.0
                break
            }

            progress <- -beta[inc_idx]/(new_beta - beta[inc_idx])
            progress <- c(progress, 1)
            
            pre_loss = objective(X2, y, new_beta, lambda, alpha)
            search_idx = order(progress)

            for(i in search_idx){
                k <- progress[i]
                
                if (k <= 0) {
                    next
                }

                if(k >= 1){
                    break
                }

                tmp_beta <- beta[inc_idx] + (new_beta - beta[inc_idx]) * k
                inner_loss <- objective(X2, y, tmp_beta, lambda, alpha)

                if (inner_loss <= pre_loss) {
                    pre_loss <- inner_loss
                    new_beta <- tmp_beta
                    line_search <- k
                }else{
                    break
                }
            }
            
            beta[inc_idx] <- new_beta
            theta[inc_idx] <- sign(new_beta)

            #  if beta encounters zero along the line search, then remove it from active set
            remove_idx = (abs(new_beta) < .Machine$double.eps)
            if(length(remove_idx) > 0){
                beta[inc_idx[remove_idx]] <- 0
                theta[inc_idx[remove_idx]] <- 0
                active_set[inc_idx[remove_idx]] <- 0
            }
        }

        exc_idx <- seq(col_num)[active_set == 0]
        grad <- XtX[exc_idx, inc_idx] %*% beta[inc_idx] - Xty[exc_idx];
        violate_idx <- (abs(grad) > lambda * alpha)

        if(length(violate_idx) == 0){
            optimality0 <- TRUE
            if(optimality1){
                break
            }
        } else {
            active_set[exc_idx[violate_idx]] <- 1
            theta[exc_idx[violate_idx]] = - sign(grad[violate_idx])
        }
        iter <- iter + 1
    }
    return(beta)
}
