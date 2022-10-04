#' Clinical Assay
#'
#' @description Simulates a clinical Assay with 2 groups (control and treatment) before and after intervention.
#'
#'
#'
#' @param n size of samples.
#' @param mean sample mean. Same for both groups before intervention (Pre-test).
#' @param sigma sample standard error.
#' @param coefvar sample coefficient of variation.
#' @param d.cohen size effect (d-Cohen). If not given, randomly generated.
#' @param dec number of decimals for observations.
#'
#' @usage
#' cassay(n, mean = 0, sigma = 1, coefvar = NULL,
#'         d.cohen = NULL, dec = 2)
#'
#' @return
#' List containing the following components:
#' \itemize{
#'
#' \item \code{Data}: a data frame containing the samples created (Columns: Group, PreTest & PostTest).
#'
#' \item \code{Model}: linear regression model.
#' }
#'
#' @examples
#' cassay(c(10,12), mean = 115, sigma = 7.5, d.cohen= 1.5)
#' cassay(24, mean = 100, sigma = 5.1)
#'
#' @export
cassay <- function(n, mean = 0, sigma = 1, coefvar = NULL, d.cohen = NULL, dec = 2){

  if(!is.null(coefvar) && coefvar != 0 && mean == 0)
    warning("Vector of means cannot be equal to zero.")

  if(length(mean) >1)
    warning("Mean must be a numeric value.")

  n <- rep(n,2)[1:2]
  sigma <- rep(sigma, 2)[1:2]
  coefvar <- rep(coefvar,2)[1:2]

  if(length(unique(n)) > 1)
    Sample <- as.factor(rep(c(0,1), times = n))
  else
    Sample <- as.factor(rep(c(0,1), each = n[1]))

  if(!is.null(coefvar)){
    if (any(mean == 0))
      stop("Parameter coefvar is given. Vector of means cannot have zeros.")
    else
      sigma <- coefvar*abs(mean)
  }

  if(is.null(d.cohen)){
    d.cohen <- stats::runif(1,-2,2)
  }

  mean_t <- mean + d.cohen*sigma # Media + efecto

  control <- pairedm(n[1], mean = mean, sigma = sigma, dec = dec)$Data
  treatment <- pairedm(n[2], mean = c(mean,mean_t), sigma = sigma, dec = dec)$Data

  Values <- rbind(control,treatment)

  d <- data.frame(cbind(Sample,Values))
  colnames(d) <- c("Group", "PreTest", "PostTest")

  l <- stats::lm(PostTest ~ Group + PreTest + Group:PreTest, data = d)
  model <- summary(l)

  return(list(Data = d, Model = model))

}



