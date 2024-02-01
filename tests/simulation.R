require(insider)

# options(error = dump_and_quit())

num_factors <- 5

# regularization for iMF
lambda <- 5
alpha <- 1 

load('~/data/insider/simulation/simulation_tensor_dataset_0.75.RData')

setwd("../results/simulation")
dataset[is.na(dataset)] <- 0

end_idx <- 2
data <- as.matrix(dataset[ ,-c(1:end_idx)]) 
confounders <- as.matrix(dataset[ ,1:end_idx])

# object <- insider(data, as.matrix(confounders), as.integer(c(1,2)), split_ratio = 0.1, global_tol = 1e-8, sub_tol = 1e-5, tuning_iter = 30)
object <- insider(data, as.matrix(confounders), split_ratio = 0.1, global_tol = 1e-8, sub_tol = 1e-5, tuning_iter = 30)
# object <- tune(object, as.integer(5), lambda = seq(1, 10, length.out = 50), alpha = 1)
object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = paste0("insider_simulation_R", num_factors, "_fitted_object.RData"))
