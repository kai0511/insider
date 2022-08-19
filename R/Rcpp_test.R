require(dplyr)
require(methods)
require(parallel)
require(Rcpp)

# load global configuration
source('./global_conf.R')
source('./imputation_functions.R')
# source('./cdplusplus.R')
sourceCpp("./cdplusplus.R")

load('~/data/zscore_matrix_selected.RData')

zscore[is.na(zscore)] <- 0
## zscore <- as(zscore, "sparseMatrix")

splitted <- ratio_splitter(zscore, return.indicator = T)
trainset <- splitted$trainset
testset <- splitted$testset
train_indicator <- (trainset != 0)
test_indicator <- (testset != 0)

names <- rownames(zscore)
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

num_factors <- 80

global_mean <- calculate_grand_mean(trainset)

r1_factor <- matrix(generate_rand_vector(num_factors * num_disease), ncol = num_factors)
r2_factor <- matrix(generate_rand_vector(num_factors * num_tissue), ncol = num_factors)
col_factor <- matrix(generate_rand_vector(num_factors * ncol(trainset)), ncol = num_factors)


trainset <- trainset - global_mean * train_indicator

loss <- 0.0; last_loss <- Inf; rmse <- 0.0; optimal_rmse = Inf
num_col <- ncol(trainset)
col_loop <- seq(num_col)


cfd1 <- as.integer(confounders[[1]])	
cfd2 <- as.integer(confounders[[2]])

r1_factor <- fitted_obj[[1]]
r2_factor <- fitted_obj[[2]]
col_factor <- fitted_obj[[3]]


or1k <- fitted_obj[[4]]
or2k <- fitted_obj[[5]]
ock <- fitted_obj[[6]]





result <- coordinate_desc_plus(trainset, testset, r1_factor, r2_factor, col_factor, cfd1, cfd2, reg_row_factor_ibias, reg_col_factor_ibias, alpha)

r1_factor <- fitted_obj[[1]]
r2_factor <- fitted_obj[[2]]
col_factor <- fitted_obj[[3]]
or1k <- fitted_obj[[4]]
or2k <- fitted_obj[[5]]
ock <- fitted_obj[[6]]
residual  <- result[[7]]
ind_matrix <- (trainset != 0)

cppFunction("arma::mat mult(arma::mat m) {
    
    change(m);
    cout << m << endl;
}

void change(arma::mat m) {
    m(0,0) = -1;
}", depends = "RcppArmadillo")



cppFunction("#include <RcppArmadillo.h>
using namespace arma;
double precision_debug(const mat& trainset, const mat& residual, 
                       const mat& r1_factor, const mat& r2_factor, const mat& col_factor, 
                       const vec& r1k, const vec& r2k, const vec& ck, 
                       const uvec& cfd1, const uvec& cfd2, const int k = 60, 
                       const double& r_lambda, const double& alpha) {
    vec pre_r1k = r1k;
    uvec zeros_idx = find(trainset == 0);
    uvec nonzeros_idx = find(trainset);
    uvec unique_cfd1 = unique(cfd1);
    update_row(zeros_idx, residual, r1k, r2k, ck, cfd1, cfd2, r_lambda, alpha);
}", depends = "RcppArmadillo")


cppFunction("Eigen::VectorXd mult(Eigen::VectorXd v) {
    Eigen::VectorXd m = (v.array() > 0).select;
    return m;
}", depends = "RcppEigen")

cppFunction("#include <iostream>
Eigen::VectorXi v = Eigen::VectorXi::Random(4);
cout << 'Here is the vector v:\n';
for(auto x : v) cout << x << ' ';
cout << '\n';", depends = "RcppEigen")



