ratio_splitter <-function(data, ratio = 0.1, rm.na.col = T){
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
