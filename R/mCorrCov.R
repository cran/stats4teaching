#' Correlation & Covariance matrices.
#' @description Given a correlation matrix and vector of standard deviations (or vector of means and vector of variation coefficients) returns a covariance matrix.
#'
#' @param mcorr a (non-empty) numeric correlation matrix.
#' @param sigma an optional vector of standard deviations.
#' @param mu    an optional vector of means.
#' @param coefvar an optional vector of coefficients of variation.
#'
#' @usage
#' mCorrCov(mcorr, sigma = 1, mu = NULL, coefvar = NULL)
#'
#' @details \code{coefvar} = \code{sigma}/\code{mu}.
#'
#'    If \code{sigma}, \code{mu} or \code{coefvar} are not specified, it´s assumed that default values for standard error's are 1. Length of standard error's is created using number of rows of correlation matrix.
#'    It's necessary to provide \code{sigma} or \code{mu} and \code{coefvar} (both) in order to obtain a desired covariance matrix.
#'
#'    Length of vectors is taken using \code{rep}. Pay attention if vectors don't have same length!
#'
#' @return \code{mCorrCov} gives the covariance matrix for a specified correlation matrix.
#'
#'
#' @examples
#' A <- matrix(c(1,2,2,1), nrow = 2, byrow = TRUE)
#' mCorrCov(A)
#'
#' B <- matrix(c(1,0.8,0.7,0.8,1,0.55,0.7,0.55,1), nrow = 3, byrow = TRUE)
#' mCorrCov(B,mu = c(2,3.5,1), coefvar = c(0.3,0.5,0.7))
#'
#' @export
mCorrCov <- function(mcorr, sigma = 1, mu = NULL, coefvar = NULL){

  k <- max(length(sigma),length(mu),length(coefvar),dim(mcorr)[1])

  if ( !is.corrmatrix(mcorr) )
    warning( "Given matrix is not a correlation matrix" )

  else if( any(sigma <= 0 ))
    warning( "Given vector of standard deviations has some element less or equal to 0")

  else if ( !is.null(coefvar) && (any(coefvar <= 0) || any(coefvar >1)))
    warning( "Coefficient of variation out of bound (0,1)" )

  else if (dim(mcorr)[1] < k)
    warning( "Dimensions of matrix of correlations are not correct" )

  else{
    mu <- rep(mu,k)[1:k]
    coefvar <- rep(coefvar,k)[1:k]

    if(!is.null(mu) && !is.null(coefvar)){
      sigma <- coefvar*abs(mu)
    }
    else{
      sigma <- rep(sigma,k)[1:k]
    }

    D <- diag(sigma)

    return(D%*%mcorr%*%D)

  }


}





