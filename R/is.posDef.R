#' Positive definited matrices
#'
#' @description Checks if a given matrix is positive definited
#' @param matrix a (non-empty) numeric matrix of data values.
#'
#' @return A logical value: True/False.
#'
#' @examples
#' A <- matrix(c(1,2,2,1), nrow = 2, byrow = TRUE)
#' is.posDef(A)
#'
#' B <- matrix(c(1,2,3,3,1,2,1,2,1), nrow = 3, byrow = TRUE)
#' is.posDef(B)
#'
#' @export
is.posDef <- function(matrix){ #comprobacion matriz es definida positiva
  if ( isSymmetric(matrix) && all(eigen(matrix)$values > 0) )
    return(TRUE)

  else
    return(FALSE)
}


