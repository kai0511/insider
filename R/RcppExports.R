# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

coordinate_descent <- function(X, y, wstart, lambda, alpha, XtX, Xty, tol = 1e-5) {
    .Call(`_insider_coordinate_descent`, X, y, wstart, lambda, alpha, XtX, Xty, tol)
}

strong_coordinate_descent <- function(X, y, wstart, lambda, alpha, XtX, Xty, tol = 1e-5) {
    .Call(`_insider_strong_coordinate_descent`, X, y, wstart, lambda, alpha, XtX, Xty, tol)
}

strong_feature_sign <- function(X, y, wstart, lambda, alpha, XtX, Xty, max_iter = 1000L) {
    .Call(`_insider_strong_feature_sign`, X, y, wstart, lambda, alpha, XtX, Xty, max_iter)
}

optimize <- function(data, cfd_factors, column_factor, cfd_indicators, train_indicator, latent_dim, lambda = 1.0, alpha = 0.1, tuning = 1L, global_tol = 1e-10, sub_tol = 1e-5, max_iter = 10000L) {
    .Call(`_insider_optimize`, data, cfd_factors, column_factor, cfd_indicators, train_indicator, latent_dim, lambda, alpha, tuning, global_tol, sub_tol, max_iter)
}

