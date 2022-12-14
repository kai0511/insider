% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insider.R
\name{tune}
\alias{tune}
\title{tune hyperparameters for the insider model and save tuning results in the local}
\usage{
tune(object, latent_dimension = NULL, lambda = 1, alpha = 0.1)
}
\arguments{
\item{object}{an insider object}

\item{latent_dimension}{a vector of ranks for selection. When tuning the rank, the default lambda and alpha are 1 and 0.1, respectively.}

\item{lambda}{a vector of l2 penalty}

\item{alpha}{a vector of l1 penalty, ranging from 0.1 to 1.}
}
\value{
tuning results
}
\description{
tune hyperparameters for the insider model and save tuning results in the local
}
\examples{
object <- tune(object, latent_dimension = as.integer(num_factors), lambda = seq(1, 50, by = 5), alpha = seq(0.1, 0.6, by = 0.1))

}
