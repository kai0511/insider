
#' dump objects when encounter errors for easy debugging
#' @export
dump_and_quit <- function() {
  # Save debugging info to file last.dump.rda
  dump.frames(to.file = TRUE)
  # Quit R with error status
  q(status = 1)
}

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

init_parameters <- function(size, init_mean = 0.0, init_std = 0.001) {
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

#' matrix splitting via element-wise sampling without replacement
#'
#' @param data matrix for splitting 
#' @param ratio a proportion of elements from data considered as testset
#' @param rm.na.col whether remove columns with all zeros
#' 
#' @return a list 
#' @export
#'
#' @examples

ratio_splitter <- function(data, ratio = 0.1, rm.na.col = T){
    # default ratio for test is 0.1, that is, 10% obs. will be randomly assigned to testset
    
    train_indicator <- matrix(T, nrow = nrow(data), ncol = ncol(data))
    testset <- matrix(0, nrow = nrow(data), ncol = ncol(data))

    test_idx <- sample(seq(length(data)), floor(length(data) * ratio), replace = F)
    testset[test_idx] <- data[test_idx]
    data[test_idx] <- 0

    train_indicator[test_idx] <- F

    num_per_col <- apply(data, 2, function(x) sum(x != 0))
    cat(paste0('number of all zero columns removed: ',  sum(num_per_col == 0)), '\n')
    if (rm.na.col) {
        return(list(trainset = data[,num_per_col != 0], 
                    testset = testset[, num_per_col != 0],
                    train_indicator = train_indicator[,num_per_col != 0]))
    } else{
        return(list(trainset = data, 
                    testset = testset,
                    train_indicator = train_indicator))
    }
}

is_converged <- function(loss, last_loss, iter, learner, thres = 1e-8, verbose = T){
    delta_loss <- last_loss - loss

    if (verbose) {
        cat(learner, " iter ", iter,  ": loss = ", loss, ", delta_loss = ", delta_loss, '\n')
    }

    if(is.na(loss) | is.infinite(loss)){
        cat("Loss = NaN or Infinity: current settings does not fit! Change the settings and try again!")
    }
    return (abs(delta_loss)/loss < thres)
}
