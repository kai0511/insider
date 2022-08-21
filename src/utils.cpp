// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


#include "../inst/include/insider_types.h"
#include <iostream>
#include "utils.h"

using std::pow;

double objective(const mat& X, const vec& y, const vec& beta,
                 const double& lambda, const double& alpha){
    vec residual = y - X * beta;
    double loss = sum(square(residual))/2;
    double l2l1_reg = (1 - alpha) * lambda * sum(square(beta))/2 + alpha * lambda * sum(abs(beta));
    loss += l2l1_reg;
    return loss;
}

double compute_loss(const vec& residual, const vec& beta, const double& lambda, const double& alpha){
    double loss = sum(square(residual))/2 + (1 - alpha) * lambda * sum(square(beta))/2 + alpha * lambda * sum(abs(beta));
    return loss;
}


void predict(const mat& row_factor, const mat& column_factor, mat& predictions) {
    predictions = row_factor * column_factor;
}

void evaluate(const mat& data, const mat& predictions, mat& residual, const uvec& train_idx, const uvec& test_idx, 
              double& sum_residual, double& train_rmse, double& test_rmse, const int& tuning, 
              const int& iter = 0, const int& verbose = 1){

    residual = data - predictions;
    sum_residual = sum(square(residual.elem(train_idx)));
    train_rmse = std::sqrt(sum_residual/train_idx.n_elem);

    if(tuning == 1){
        test_rmse = std::sqrt(mean(square(residual.elem(test_idx))));
    }  

    if (verbose == 1){
        cout << "ibiasedMF with L1_penalty iter " << iter << ": train rmse = " << train_rmse << '\n' << endl;

        if(tuning == 1){
            cout << "ibiasedMF with L1_penalty iter " << iter << ": train rmse = " << test_rmse << '\n' << endl;
        }
    }
}

double compute_loss(const mat& cfd_factor, const mat& column_factor, const double& lambda, const double& alpha, 
                    double& sum_residual, const int& verbose = 1){

    // l2 penalty 
    double row_reg = lambda * pow(norm(cfd_factor, "F"), 2);
    double col_reg = lambda * (1 - alpha) * pow(norm(column_factor, "F"), 2);
    
    //  l1 penalty
    double l1_reg = lambda * alpha * sum(sum(abs(column_factor), 1));

    double loss = sum_residual/2 + row_reg/2 + col_reg/2 + l1_reg;

    if(verbose == 1) {
        cout << "total_residual" << '\t' << sum_residual/2 << ";" << '\n'
             << "row_reg_loss:" << '\t' << row_reg/2 << ";" << '\n'
             << "col_reg_loss:" << '\t' << col_reg/2 << ";" << '\n'
             << "l1_reg_loss:" << '\t' << l1_reg << "." << endl;
    }
    return loss;
}
