#' Correlation matrix
#'
#' @description Checks if a given matrix is a correlation matrix for non-degenerate distributions.
#'
#' @param matrix a (non-empty) numeric matrix of data values.
#'
#' @return A logical value: True/False.
#'
#'
#' @examples
#'
#' m1<-matrix(c(1,2,2,1),nrow = 2,byrow = TRUE)
#' is.corrmatrix(m1)
#'
#' m2<-matrix(c(1,0.8,0.8,1),nrow = 2,byrow = TRUE)
#' is.corrmatrix(m2)
#'
#' m3<-matrix(c(1,0.7,0.8,1),nrow = 2,byrow = TRUE)
#' is.corrmatrix(m3)
#'
#' @export
is.corrmatrix <- function(matrix){ #comprobacion matriz correlacion
  if( any(diag(matrix) != 1) || !is.posDef(matrix) || any(abs(matrix) > 1) )
    return(FALSE)
  else TRUE

}
