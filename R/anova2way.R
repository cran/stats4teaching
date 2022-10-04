#' Two-Way ANOVA
#'
#' @description \code{anova2way} returns multivariate data in order to compute analysis of variance with 2 factors.
#'
#'
#' @param k number of levels Factor I. By default k=2.
#' @param j number of levels Factor II. By default j=2.
#' @param n number of elements in each group (k,j).
#' @param mean vector of means.
#' @param sigma vector of standard deviations.
#' @param coefvar an optional vector of coefficients of variation.
#' @param method post-hoc method applied. There are five possible choices: ``\code{Tukey}``, ``\code{LSD}``, ``\code{Dunnett}``, ``\code{Bonferroni}``, ``\code{Scheffe}``. Can be specified just the initial letter.
#' @param conf.level confidence level of the interval.
#' @param dec number of decimals for observations.
#'
#' @usage
#' anova2way(k =2 , j = 2, n,  mean = 0, sigma = 1,
#'           coefvar = NULL, method = c("Tukey", "LSD", "Dunnett", "Bonferroni", "Scheffe"),
#'           conf.level = 0.95, dec = 2)
#'
#' @return A list containing the following components:
#' \itemize{
#'
#'  \item \code{Data}: a data frame containing the samples created.
#'
#'
#'  \item \code{Size.effect}: size effect for each factor and interaction.
#'
#'  \item \code{Significance/Test Post-Hoc}: significance for each factor and interaction and test Post-Hoc for each factor.
#'
#' }
#' @examples
#'
#' anova2way(k=3, j=2, n=c(3,4,4,5,5,3), mean = c(1,4,2.5,5,6,3.75), sigma = c(1,1.5))
#'
#' @export

anova2way <- function(k = 2, j = 2, n, mean = 0, sigma = 1, coefvar =  NULL, method = c("Tukey", "LSD", "Dunnett", "Bonferroni", "Scheffe"), conf.level = 0.95, dec = 2){
   # k = number levels of Factor I
  # j = number levels of Factor II
  # n = sample size for observations.
  if(length(k) > 1)
    stop("Number of elements in each group must be equal")
  if((length(mean) || length(sigma)) > k*j)
    stop("Length of vectors are not compatible with levels of factors")
  q <- max(length(mean), length(sigma), length(coefvar), k*j,n)
  n <- rep(n,q)[1:q]
  mean <- rep(mean, q)[1:q]
  coefvar <- rep(coefvar, q)[1:q]
  if(!is.null(coefvar))
    sigma <- coefvar * abs(mean)
  else
    sigma <- rep(sigma, q)[1:q]
  # FactorI
  FactorI <- c()
  i = 1
  for(m in 1:k){
    b <- NULL
    for(l in 1:j){
      a <- rep(m, times = n[i])
      FactorI <- c(FactorI, a)
      i = i+1
      b = i
    }
  }
  FactorI <- as.factor(FactorI)
  # FactorII
  FactorII <- c()
  i = 1
  for(m in 1:k*j){
    b <- NULL
    for(l in 1:j){
      a <- rep(l, times = n[i])
      FactorII <- c(FactorII, a)
      i = i+1
      b = i
    }
  }
  FactorII <- as.factor(FactorII)
  #Generating data.
  Values <- c()
  for(i in 1:(k*j)){
    gen <- generator(n[i], mean[i], sigma[i], dec = dec)
    Values <- c(Values, gen$Sample)
  }

  d <- data.frame(FactorI, FactorII, Values)
  #ANOVA
  anova <- stats::aov(Values ~ FactorI*FactorII, data = d)
  anova <- car::Anova(anova, type = "II")
  p.value <- round(anova$`Pr(>F)`[1:3], digits = 4)
  p.valueS <- c()
  for(i in 1:3){
    if(p.value[i] <= 0.001)
      p.valueS[i] <- "< 0.001"
    else
      p.valueS[i] <- round(p.value[i], digits = 3)
  }
  alpha <- 1 - conf.level
  # Levene's test for variances homogeneity. Informative.
  levene_H <- car::leveneTest(Values ~ FactorI*FactorII, data = d)
  # Printing ANOVA Table (homogeneous decimals)
  dig <- options(digits = 4)
  on.exit(options(dig))
  printA <- knitr::kable(anova, digits = 4, format = "simple")
  eta <- round(rstatix::eta_squared(anova),4) # Size effect.

  #Post-Hoc test
  if(is.null(method) || method[1] == "Tukey" || method[1] == "T"){
    postHocI <- asbio::tukeyCI(Values, FactorI, conf.level)
    postHocII <- asbio::tukeyCI(Values, FactorII, conf.level)
  }
  else if(method[1] == "LSD" || method[1] == "L"){
    postHocI <- asbio::lsdCI(Values, FactorI, conf.level)
    postHocII <- asbio::lsdCI(Values, FactorII, conf.level)
  }
  else if(method[1] == "Dunnett" || method[1] == "D"){
    postHocI <- asbio::dunnettCI(Values, FactorI, conf.level)
    postHocII <- asbio::dunnettCI(Values, FactorII, conf.level)
  }
  else if(method[1] == "Bonferroni" || method[1] == "B"){
    postHocI <- asbio::bonfCI(Values, FactorI, conf.level)
    postHocII <- asbio::bonfCI(Values, FactorII, conf.level)
  }
  else{
    postHocI <- asbio::scheffeCI(Values, FactorI, conf.level)
    postHocII <- asbio::scheffeCI(Values, FactorII, conf.level)
  }

  if(p.value[1] < alpha){
    SignificanceI <- paste(p.valueS[1],", significative ANOVA for Factor I", sep = "")
  }
  else{
    SignificanceI <- paste(p.valueS[1],", not significative ANOVA for Factor I", sep = "")
  }
  if(p.value[2] < alpha){
    SignificanceII <- paste(p.valueS[2],", significative ANOVA for Factor II", sep = "")
  }
  else{
    SignificanceII <- paste(p.valueS[2],", not significative ANOVA for Factor II", sep = "")
  }
  if(p.value[3] < alpha){
    SignificanceIII <- paste(p.valueS[3],", significative ANOVA for interaction", sep = "")
  }
  else
    SignificanceIII <- paste(p.valueS[3],", not significative ANOVA for interaction", sep = "")
  return(list(Data = d, leveneTest = levene_H, Anova = printA, Size.effect = eta,
              FactorI = list(Significance = SignificanceI, PostHoc.Test = postHocI),
              FactorII = list(Significance = SignificanceII, PostHoc.Test = postHocII),
              Interaction = SignificanceIII))
}
