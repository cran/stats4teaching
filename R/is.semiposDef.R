#' Semi-Positive definited matrices
#'
#' @description Checks if a given matrix is semi-positive definited.
#' @param matrix a (non-empty) numeric matrix of data values.
#'
#' @return A logical value: True/False.
#'
#' @examples
#' A<-matrix(c(2.2,1,1,3), nrow = 2, byrow = TRUE)
#' is.semiposDef(A)
#'
#' B<-matrix(c(1,2,3,3,1,2,1,2,1), nrow = 3, byrow = TRUE)
#' is.semiposDef(B)
#'
#' @export
is.semiposDef <- function(matrix){
  if ( isSymmetric(matrix) && all(eigen(matrix)$values >= 0))
    return(TRUE)

  else
    return(FALSE)
}
