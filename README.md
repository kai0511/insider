# INSIDER

INSIDER is a computational tool for interpretable sparse matrix decomposition for bulk RNA expression data analysis. It decomposes variation from multiple biological conditions into a shared low-rank latent space. Our approach enables various downstream analysis. See the preprint in References section for details of downstream analysis.

## Dependencies
This package relies on the following R packages. 
```{r}
Rcpp (>= 1.0.7), RcppArmadillo, dplyr, MASS
```

## Installation
INSIDER uses **openmp** to support parallel computing. To enjoy the benefit, one can install **openmp** before installation. This package can be installed via the following ways.

### Install via devtools
```{r}
install.packages("devtools")
install_github("kai0511/insider")
```
### Install locally
Download the zip file from GitHub, unzip it, and go into the directory. Then, it can be installed by running the following commands
```{Shell}
R CMD build .
R CMD INSTALL insider_1.0.tar.gz 
```

## Usage

* Data preparation
```{r}
require(insider)

load('~/data/multi_dimensional_datasets/ageing_dataset_annotated_with_phenotypes_filtered.RData')
dataset <- dataset[,-1]

end_idx <- 3  # The end index for covariate matrix
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)  # log transformed expression data matrix
data[is.na(data)] <- 0 # cast NAs to zeros

confounders <- as.matrix(dataset[ ,1:end_idx])   # matrix for covariates ()
```
* Create INSIDER object
```{r}
object <- insider(data, confounders, interaction_idx = as.integer(c(1,2)), split_ratio = 0.1, global_tol = 1e-9, sub_tol = 1e-5, tuning_iter = 30)
```
It needs the following arguments:
1. *data*: A log-transformed expression data matrix;
2. *confounder*: A confounder matrix;
3. *interaction_idx*: An integer vector for indices of confounders to induce interaction. For example, as.integer(c(1,2)) means to consider interaction between the first and second variables;
4. *split_ratio*: define the proportion of elements in the data matrix used as test set for hyperparameter tuning.  
5. *global_tol*: defines convergence tolerance for INSIDER. Note INSIDER check convergence every 10 iterations;
6. *sub_tol*: defines the convergence criteria for elastic net problems;
7. *tuning_iter*: the number of iterations to run for each try of hyperparameter combinations.
8. *max_iter*: the maximum number of iterations. When it is reached, iteration will terminate even if the global convergence criteria do not meet.

* Tune hyperparameters
```{r}
object <- tune(object, latent_dimension = as.integer(seq(10, 30, by = 2)), lambda = seq(1, 20, by = 2), alpha = c(0.2, 0.3, 0.4, 0.5))
```
It needs the following arguments:
1. *object*: An INSIDER object created with the above arguments;
2. *latent_dimension*: A integer vector from which the rank of latent dimension is chosen;
3. *lambda*: A numeric vector from which the tuning parameter *lambda* is selected;
4. *alpha*: A numeric vector from which the tuning parameter *alpha* is selected;

* After parameter tuning, the results for tuning will be saved in the current directory. One chose the combination of hyperparameters with the lowest RMSE on test, and fit INSIDER with it.
```{r}
# selected hyperparameters for INSIDER
num_factors <- 23
lambda <- 10
alpha <- 0.4
object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = paste0("insider_ageing_R", num_factors, "_fitted_object.RData"))
```