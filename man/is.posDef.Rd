% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.posDef.R
\name{is.posDef}
\alias{is.posDef}
\title{Positive definited matrices}
\usage{
is.posDef(matrix)
}
\arguments{
\item{matrix}{a (non-empty) numeric matrix of data values.}
}
\value{
A logical value: True/False.
}
\description{
Checks if a given matrix is positive definited
}
\examples{
A <- matrix(c(1,2,2,1), nrow = 2, byrow = TRUE)
is.posDef(A)

B <- matrix(c(1,2,3,3,1,2,1,2,1), nrow = 3, byrow = TRUE)
is.posDef(B)

}
