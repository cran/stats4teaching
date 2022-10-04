#' One-Way ANOVA
#'
#' @description \code{anova1way} is used to generate multivariate data in order to compute analysis of variance with 1 factor. It provides balanced and unbalanced ANOVA (as long as homogeneity of variances is satisfied. In other case it is provided Welch test).
#'
#'
#'
#' @param k number of levels. By default k = 3.
#' @param n size of samples.
#' @param mean vector of means.
#' @param sigma vector of standard deviations.
#' @param coefvar an optional vector of coefficients of variation.
#' @param method post-hoc method applied. There are five possible choices: "\code{Tukey}", "\code{LSD}", "\code{Dunnett}", "\code{Bonferroni}", "\code{Scheffe}". Can be specified just the initial letter.
#' @param conf.level confidence level of the interval.
#' @param dec number of decimals for observations.
#'
#' @usage
#' anova1way(k = 3,n , mean = 0, sigma = 1,
#'           coefvar = NULL, method = c("Tukey", "LSD", "Dunnett", "Bonferroni", "Scheffe"),
#'           conf.level = 0.95, dec = 2)
#'
#'
#' @details If \code{mean} or \code{sigma} are not specified it is assumed the default values of \code{0} and \code{1}.
#'
#'     If \code{coefvar} (= \code{sigma}/\code{mean}) is specified, function omits \code{sigma}.
#'
#'     Number of samples is choosen by \code{k} (by default k = 3). Therefore, if the others parameters (\code{n}, \code{mean}, \code{sigma}, \code{coefvar}) have not same length, function \code{rep} will be used. Pay attention if vectors dont have same length.
#'
#'     Moreover, not only gives samples for each level, but also the ANOVA table and post-hoc test (in case of significance). By default \code{conf.level} = 0.95 and Tukey method is used. If the homogeneity of variances is not verified (using Bartlett test), the Welch test is performed.
#'
#' @return List containing the following components:
#'
#' \itemize{
#'
#'  \item \code{Data}: a data frame containing the samples created.
#'
#'  \item \code{Anova}: anova fitted model.
#'
#'  \item \code{Significance}: significance of the factor.
#'
#'  \item \code{Size.effect}: size effect of the factor.
#'
#'  \item \code{Test Post-Hoc}: test Post-Hoc.
#' }
#'
#' @examples
#' anova1way(k=4,n=c(40,31,50),mean=c(55,52,48,59),coefvar=c(0.12,0.15,0.13),conf.level = 0.99)
#'
#' anova1way(k=3,n=15,mean=c(10,15,20),sigma =c(1,1.25,1.1),method ="B")
#'
#'
#' @export

anova1way <- function(k = 3, n ,mean = 0 , sigma = 1, coefvar=  NULL, method = c("Tukey", "LSD", "Dunnett", "Bonferroni", "Scheffe"), conf.level = 0.95, dec = 2){
  # Balanced and not balanced ANOVA


  q <- max( length(n), length(mean) , length(sigma) , length(coefvar) )

  if(k < q) #Some parameter is longer than number of groups.
    stop("Number of groups and parameters incompatible")

  n <- rep(n, k)[c(1:k)]
  mean <- rep(mean, k)[c(1:k)]
  coefvar <- rep(coefvar, k)[1:k]

  if( is.null(coefvar) )
    sigma <- rep(sigma, k)[c(1:k)]
  else
    sigma <- coefvar * abs(mean)
  # Same length for all parameters.

  gen <- generator(n, mean, sigma, dec = dec)

  Values <- c()
  if(length(unique(n)) > 1){ # Sample sizes are different.
    for(i in 1:k)
      Values <- c(Values, gen[[i]]) # List format.
  }

  else{ # Sample sizes are equal.
    for(i in 1:k)
      Values <- c(Values, gen$Samples[,i]) # data.frame format.
  }


  Sample <- as.factor(rep(c(1:k), times = n))
  samples <- data.frame(Values, Sample)

  alpha <- 1 - conf.level

  #Bartlett test for homogeneity of variances.
  bartlett_H <- stats::bartlett.test(Values ~ Sample, data = samples)
  p.value_bartlett <- bartlett_H$p.value

  if(p.value_bartlett  < alpha){
    # Welch ANOVA
    w_test <- stats::oneway.test(Values ~ Sample, data = samples, var.equal = FALSE)

    if(w_test$p.value <= 0.001)
      p.valueW <- "< 0.001"
    else
      p.valueW <- round(w_test$p.value, digits = 4)

    if(is.null(method) || method[1] == "Tukey" || method[1] == "T")
      postHoc <- asbio::tukeyCI(Values,Sample, conf.level)

    else if(method[1] == "LSD" || method[1] == "L")
      postHoc <- asbio::lsdCI(Values,Sample, conf.level)

    else if(method[1] == "Dunnett" || method[1] == "D")
      postHoc <- asbio::dunnettCI(Values,Sample, conf.level)

    else if(method[1] == "Bonferroni" || method[1] == "B")
      postHoc <- asbio::bonfCI(Values,Sample, conf.level)

    else
      postHoc <- asbio::scheffeCI(Values,Sample, conf.level)

    if(w_test$p.value < alpha) {
      return(list(Data = samples, Anova = "Cannot assume homogeneity of variances (by  Bartlett test)", Welch.ANOVA = w_test, Significance = paste(p.valueW,", significative Welch ANOVA", sep = ""), PostHoc.Test = postHoc))
    }

    else{

      return(list(Data = samples, Anova = "Cannot assume homogeneity of variances (by  Bartlett test)", Welch.ANOVA = w_test, Significance = paste(p.valueW,", not significative Welch ANOVA", sep = ""), PostHoc.Test = postHoc))
    }

  }


  else{

    # Balanced or unbalanced designs
    anova <- stats::aov(Values ~ Sample, data = samples)
    s <- summary(anova)

    p.value <- s[[1]]$`Pr(>F)`[1]

    if(p.value <= 0.001)
      p.valueS <- "< 0.001"
    else
      p.valueS <- round(s[[1]]$`Pr(>F)`[1], digits = 4)

    dig <- options(digits = 4)
    on.exit(options(dig)) # Printing ANOVA Table (homogeneous decimals)
    printAnova <- knitr::kable(s[[1]], digits = 4, format = "simple")


    eta <- round(rstatix::eta_squared(anova)[[1]], 4) #Size effect

    if(is.null(method) || method[1] == "Tukey" || method[1] == "T")
      postHoc <- asbio::tukeyCI(Values,Sample, conf.level)

    else if(method[1] == "LSD" || method[1] == "L")
      postHoc <- asbio::lsdCI(Values,Sample, conf.level)

    else if(method[1] == "Dunnett" || method[1] == "D")
      postHoc <- asbio::dunnettCI(Values,Sample, conf.level)

    else if(method[1] == "Bonferroni" || method[1] == "B")
      postHoc <- asbio::bonfCI(Values,Sample, conf.level)

    else
      postHoc <- asbio::scheffeCI(Values,Sample, conf.level)

    if(p.value < alpha){
      return(list(Data = samples, Anova = printAnova, Significance = paste(p.valueS,", significative ANOVA", sep = ""), Size.effect = eta, PostHoc.Test = postHoc))
    }

    else{
      return(list(Data = samples, Anova = printAnova, Significance = paste(p.valueS,", not significative ANOVA", sep = ""), PostHoc.Test = postHoc))
    }
  }


}
