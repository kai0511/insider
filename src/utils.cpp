// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


#include "../inst/include/insider_types.h"
#include <iostream>
#include "utils.h"

using std::pow;

template <typename T>
inline bool rows_equal(const T& lhs, const T& rhs, double tol = 1e-8) {
    return approx_equal(lhs, rhs, "absdiff", tol);
}

mat unique_rows(const mat& m) {
    uvec flag = zeros<uvec>(m.n_rows);
    for (uword i = 0; i < m.n_rows; i++) {
        for (uword j = i + 1; j < m.n_rows; j++) {
            if (rows_equal(m.row(i), m.row(j))) { flag(j) = 1; break; }
        }
    }
    return m.rows(find(flag == 0));
}

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

void evaluate(mat& residual, const uvec& train_idx, const uvec& test_idx, 
              double& sum_residual, double& train_rmse, double& test_rmse, const int& tuning, 
              const int& iter = 0, const int& verbose = 1){

    // residual = data - predictions;
    if(tuning == 0){
        sum_residual = accu(sum(square(residual), 1));
        train_rmse = std::sqrt(sum_residual/(residual.n_cols * residual.n_rows));
    }else{
        sum_residual = sum(square(residual.elem(train_idx)));
        train_rmse = std::sqrt(sum_residual/train_idx.n_elem);
        test_rmse = std::sqrt(mean(square(residual.elem(test_idx))));
    } 

    if (verbose == 1){
        cout << "insider iter " << iter << ": train rmse = " << train_rmse << endl;

        if(tuning == 1){
            cout << "insider iter " << iter << ": test rmse = " << test_rmse << endl;
        }
    }
}

double compute_loss(const field<mat>& cfd_factor, const mat& column_factor, const double& lambda, const double& alpha, 
                    double& sum_residual, const int& verbose = 1){

    // l2 penalty 
    double row_reg = 0.0;
    for(unsigned int i = 0; i < cfd_factor.n_elem; i++){
        row_reg += lambda * pow(norm(cfd_factor(i), "F"), 2);
    }
    
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
