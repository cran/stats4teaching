% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.covmatrix.R
\name{is.covmatrix}
\alias{is.covmatrix}
\title{Covariance matrix}
\usage{
is.covmatrix(matrix)
}
\arguments{
\item{matrix}{a (non-empty) numeric matrix of data values.}
}
\value{
A logical value: True/False.
}
\description{
Checks if a given matrix is a covariance matrix for non-degenerate distributions.
}
\examples{

m1 <- matrix(c(2,1.5,1.5,1), nrow = 2, byrow = TRUE)
is.covmatrix(m1)

m2 <- matrix(c(1,0.8,0.8,1), nrow = 2, byrow = TRUE)
is.covmatrix(m2)

m3 <- matrix(c(1,0.7,0.8,1), nrow = 2, byrow = TRUE)
is.covmatrix(m3)

}
