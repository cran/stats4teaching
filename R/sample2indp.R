#' Independent normal data
#'
#' @description Generates two normal independent samples. It also provides Cohen's effect and T-Test.
#'
#'
#' @param n vector of size of samples.
#' @param mean vector of means.
#' @param sigma vector of standard deviations.
#' @param coefvar an optional vector of coefficients of variation.
#' @param alternative a character string specifying the alternative hypothesis for T-Test. meanst be one of ``two.sided`` (default), ``greater`` or ``less``. Can be specified just the initial letter.
#' @param delta true value of the difference in means.
#' @param conf.level confidence level of the interval. It determines level of significance for comparing variances.
#' @param dec number of decimals for observations.
#'
#' @usage
#' sample2indp(n , mean = 0, sigma = 1, coefvar = NULL,
#'             alternative = c("two.sided", "less", "greater"), delta = 0,
#'             conf.level = 0.95, dec = 2)
#'
#' @details If \code{mean} or \code{sigma} are not specified it's assumed the default values of \code{0} and \code{1}.
#'
#'    \code{n} is a vector, so it's possible to generate samples with same or different sizes.
#'
#'    If \code{coefvar} is given, \code{sigma} is omitted. Vector of means cannot have any 0.
#'
#' @return A list containing the following components:
#'
#' \itemize{
#'
#'  \item \code{Data}: a data frame containing the samples created.
#'
#'
#'  \item \code{T.Test}: a t-test of the samples.
#'
#'  \item \code{Power}: power of the test.
#'
#' }
#'
#' @examples
#' sample2indp(c(10,12),mean = c(2,3),coefvar = c(0.3,0.5), alternative = "less", delta = -1)
#'
#' sample2indp(8,sigma = c(1,1.5), dec = 3)
#'
#' @export

sample2indp <- function(n, mean = 0, sigma = 1, coefvar = NULL, alternative = c("two.sided", "less", "greater"), delta = 0, conf.level = 0.95, dec = 2){

  if(any(length(n) > 2, length(unique(n)) > 2,length(mean) > 2, length(sigma) > 2, length(coefvar) > 2))
    stop("More of 2 samples are being created")

  else if(!is.null(sigma) && !is.vector(sigma))
    stop("Sigma meanst be a vector")

  k <- 2

  mean <- rep(mean,k)[c(1:k)]
  coefvar <- rep(coefvar,k)[c(1:k)]

  if(!is.null(coefvar)){
    if (any(mean == 0))
      stop("Parameter coefvar is given. Vector of means cannot have zeros.")
    else
      sigma <- coefvar*abs(mean)
  }


  else
    sigma <- rep(sigma,k)[c(1:k)]

  n <- rep(n,k)[c(1:k)]  #Parameters with same lenght.

  gen <- generator(n, mean, sigma, dec = dec)

  if(length(unique(n)) == 1){ #Size of samples is equal.
    Values <- c(gen$Samples[, 1], gen$Samples[, 2]) #data.frame (Same sizes)
    Sample <- as.factor(rep(c(1,2), each = n[1]))

  }
  else{
    Values <- c(gen$X1, gen$X2) #List (Different sizes)
    Sample <- as.factor(rep(c(1,2), times = n))

  }

  samples <- data.frame(Values, Sample)


  # Variance Test

  vtest <- stats::var.test(Values ~ Sample, data = samples)
  sp <- sqrt(((n[1] - 1) * sigma[1]^2 + (n[2] - 1) * sigma[2]^2) / (n[1] + n[2] - 2)) # Pooled standard error
  d.cohen <- round((mean[1] - mean[2]) / sp, 4) #Cohen's effect

  if(abs(d.cohen) <= 0.2)
    e.cohen <- paste(d.cohen, "(small)")
  else if(abs(d.cohen) <= 0.5)
    e.cohen <- paste(d.cohen, "(moderate)")
  else
    e.cohen <- paste(d.cohen, "(big)")

  #T-Test
  alpha <- 1 - conf.level

  if(vtest$p.value >=  alpha)
    var <- TRUE
  else
    var <- FALSE

  t <- stats::t.test(Values ~ Sample, data = samples, mu = delta, alternative = alternative[1], var.equal = var, conf.level = conf.level)
  t$data.name <- c("Sample1 ","Sample2")

  # Power test

  if(length(unique(n)) == 1)
    pow <- pwr::pwr.t.test(n[1], d = d.cohen, sig.level = 1- conf.level, power = NULL, type = "two.sample", alternative = alternative[1])$power

  else
    pow <- pwr::pwr.t2n.test(n1 = n[1], n2 = n[2], d = d.cohen, sig.level = 1- conf.level, power = NULL, alternative = alternative[1])$power


  return(list(Data = samples, Cohen.effect = e.cohen, T.test = t, power = pow))
}

