#' Covariance matrix
#'
#' @description Checks if a given matrix is a covariance matrix for non-degenerate distributions.
#'
#' @param matrix a (non-empty) numeric matrix of data values.
#'
#' @return A logical value: True/False.
#'
#'
#' @examples
#'
#' m1 <- matrix(c(2,1.5,1.5,1), nrow = 2, byrow = TRUE)
#' is.covmatrix(m1)
#'
#' m2 <- matrix(c(1,0.8,0.8,1), nrow = 2, byrow = TRUE)
#' is.covmatrix(m2)
#'
#' m3 <- matrix(c(1,0.7,0.8,1), nrow = 2, byrow = TRUE)
#' is.covmatrix(m3)
#'
#' @export
is.covmatrix<-function(matrix){ #comprobacion matriz correlacion
  if( !is.posDef(matrix) || !isSymmetric(matrix) )
    return(FALSE)
  else TRUE

}

#Teoria: http://halweb.uc3m.es/esp/Personal/personas/agrane/esp/cooperacion/proyecto_mozambique_archivos/alumnos/Tema3_NormalidadMultivariante_reducido.pdf
#Teoria: https://encyclopediaofmath.org/wiki/Covariance_matrix
