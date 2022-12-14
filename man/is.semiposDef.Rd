% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.semiposDef.R
\name{is.semiposDef}
\alias{is.semiposDef}
\title{Semi-Positive definited matrices}
\usage{
is.semiposDef(matrix)
}
\arguments{
\item{matrix}{a (non-empty) numeric matrix of data values.}
}
\value{
A logical value: True/False.
}
\description{
Checks if a given matrix is semi-positive definited.
}
\examples{
A<-matrix(c(2.2,1,1,3), nrow = 2, byrow = TRUE)
is.semiposDef(A)

B<-matrix(c(1,2,3,3,1,2,1,2,1), nrow = 3, byrow = TRUE)
is.semiposDef(B)

}
