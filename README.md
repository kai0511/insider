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
Here we use a small proportion of the ageing dataset (377*44477) as data matrix for a toy example. For demonstration purpose, the dimension of example data matrix is 377 by 5000. One can download the data for analysis with INSIDER.
```{r}
require(insider)

load('ageing_example_data.RData')   # load the data and the exact path for the example data matrix depends.

head(dataset[,1:5])
  rnaseq_profile_id pid sid did X499304660
1         496100536   1   4  71  0.0000000
2         496100556   1   6  71  0.0000000
3         496100513   1   6  69  0.0000000
4         496100514   1   2  69  0.0000000
5         496100520   1   4  69  0.2626623
6         496100524   1   5  69  0.0000000

dataset <- dataset[,-1]   # The first column is rna-seq_profile_id

end_idx <- 3  # The end index for covariate matrix
data[is.na(data)] <- 0 # cast NAs to zeros
data <- log2(as.matrix(dataset[ ,-c(1:end_idx)]) + 1)  # log transformed expression data matrix

# In the example data, there are three biological variables: pid (phenotype id), sid (brain structure id), and did (donor id).
confounders <- as.matrix(dataset[ ,1:end_idx])   # matrix for biological variables
```

* Create INSIDER object
```{r}
object <- insider(data, confounders, interaction_idx = as.integer(c(1,2)), split_ratio = 0.1, global_tol = 1e-9, sub_tol = 1e-5, tuning_iter = 30)
```
It needs the following arguments:
1. *data*: A log-transformed expression data matrix;
2. *confounder*: A confounder matrix. The elements of the matrix are used as indices to extract corresponding latent representation, so its elements are integer and greater than 0;
3. *interaction_idx*: An integer vector for indices of confounders to induce interaction. For example, as.integer(c(1,2)) means to consider interaction between phenotype and brain structures, with the above example data;
4. *split_ratio*: define the proportion of elements in the data matrix used as test set for hyperparameter tuning.  
5. *global_tol*: defines convergence tolerance for INSIDER. Note INSIDER check convergence every 10 iterations, so global_tol equal 1e-9 is equivalent to the stopping criteria defined in our preprint in references.
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

> str(object)
List of 7
 $ data           : num [1:377, 1:5000] 0 0 0 0 0.336 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:5000] "X499304660" "X499304661" "X499304664" "X499304666" ...
 $ confounder     : num [1:377, 1:4] 1 1 1 1 1 1 2 2 2 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "pid" "sid" "did" "interaction_indicator"
 $ train_indicator: int [1:377, 1:5000] 1 1 1 1 0 1 1 1 1 1 ...
 $ params         :List of 4
  ..$ global_tol : num 1e-10
  ..$ sub_tol    : num 1e-05
  ..$ tuning_iter: num 30
  ..$ max_iter   : num 50000
 $ cfd_matrices   :List of 4
  ..$ factor0: num [1:2, 1:23] -0.155 0.128 0.144 -0.217 0.862 ...
  ..$ factor1: num [1:8, 1:23] -0.0301 -0.057 -0.2099 0.1279 0.0856 ...
  ..$ factor2: num [1:107, 1:23] -0.377 -0.778 -0.11 -0.552 0.251 ...
  ..$ factor3: num [1:16, 1:23] -0.000318 0.196197 -0.040399 0.031864 0.053713 ...
 $ column_factor  : num [1:23, 1:5000] 0 0.00386 0 -0.00831 0 ...
 - attr(*, "class")= chr "insider"
```

The fitted object obtained from the above command is an R list object, containing the following elements:
1. log-transformed expression data matrix;
2. confounder matrix
3. train_indicator: an indicator matrix for elements to be concluded as train set.
4. params: parameter setting for INSIDER
6. cfd_matrices: a list of low-rank representations for biological variables and interaction. One can access the low-rank representation for a specific biological variable with the index of the variable in the confounder matrix. The low-rank representation for the interaction is the last matrix of the cfd_matrices List.
7. column_factor: gene latent representation matrix of K * M, where K is the num_factors and M is the number of genes.

For downstream analysis with results from INSIDER, please refer to preprint in references. 

## References
Zhao, Kai, et al. "INSIDER: Interpretable Sparse Matrix Decomposition for Bulk RNA Expression Data Analysis." bioRxiv (2022): 2022-11.
