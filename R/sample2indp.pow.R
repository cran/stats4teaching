#' Independent normal data
#'
#' @description Generates two normal independent samples with desired power and cohen's effect.
#'
#' @param n1 first sample size.
#' @param mean vector of sample means.
#' @param s1 standard deviation for first sample.
#' @param d.cohen Cohen's effect.
#' @param power power of the test.
#' @param alternative a character string specifying the alternative hypothesis for T-Test. Must be one of ``two.sided`` (default), ``greater`` or ``less``. Can be specified just the initial letter.
#' @param delta true value of the difference in means.
#' @param conf.level confidence level of the interval.
#' @param dec number of decimals for observations.
#'
#' @usage sample2indp.pow(n1, mean = 0, s1= 1, d.cohen, power,
#'    alternative = c("two.sided", "less", "greater"), delta = 1,
#'    conf.level = 0.95, dec = 2)
#'
#' @details
#'
#'    Pooled standard deviation= \code{sp} = sqrt{((n1 - 1) sigma1^2 +(n2 - 1) sigma2^2) / (n1 + n2 - 2)}
#'
#'     \code{d.cohen} = |mean1 - mean2| / sqrt(sp)
#'
#' @return A list containing the following components:
#'
#' \itemize{
#'
#'  \item \code{Data}: a data frame containing the samples created.
#'
#'
#'  \item \code{Size}: size of each sample.
#'
#'  \item \code{T.test}: a t-test of the samples.
#'
#' }
#'
#'
#' @examples
#' sample2indp.pow(n1 = 30, mean = c(2,3), s1= 0.5, d.cohen = 0.8, power = 0.85, delta = 1)
#' sample2indp.pow(n1 = 50, mean = c(15.5,16), s1=2 , d.cohen = 0.3, power = 0.33, delta = 0.5)
#'
#' @export

sample2indp.pow <- function(n1, mean = 0, s1= 1, d.cohen, power, alternative = c("two.sided", "less", "greater"), delta = 1, conf.level = 0.95, dec = 2){

  n2 <- ceiling(pwr::pwr.t2n.test(n1 = n1, n2 = NULL, d = d.cohen , sig.level = 1 - conf.level, power, alternative = alternative[1])$n2)
  n <- c(n1,n2)

  sp <- abs((mean[2]-mean[1]))/d.cohen

  if(sp^2 * (n1 + n2 - 2) <= (n1 - 1) * s1^2)
    stop("Standard deviation of first sample and pooled standard deviation are not compatible")

  else{

    s2 <- sqrt((sp^2 * (n1 + n2 - 2) - (n1 - 1) * s1^2) / (n2 - 1))
    sigma <- c(s1, s2)

    gen <- generator(n, mean, sigma, dec = dec)

    if(length(unique(n)) == 1) #Size of samples is equal.
      Values <- c(gen$Samples[, 1], gen$Samples[, 2]) #data.frame (Same sizes)

    else
      Values <- c(gen$X1, gen$X2) #List (Different sizes)

  }


  Sample <- as.factor(rep(c(1,2), times = n))
  samples <- data.frame(Values, Sample)


  vtest <- stats::var.test(Values ~ Sample, data = samples)

  #T-Test
  alpha <- 1 - conf.level

  if(vtest$p.value >=  alpha)
    var <- TRUE
  else
    var <- FALSE

  t <- stats::t.test(Values ~ Sample, data = samples, mu = delta, alternative = alternative[1], var.equal = var, conf.level = conf.level)
  t$data.name <- c("Sample1 ","Sample2")


  return(list(Data = samples, Size = n, T.test = t))

}
