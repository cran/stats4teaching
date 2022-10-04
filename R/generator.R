#' Generation of multivariate normal data.
#'
#' @description This function generates univariate and multivariate normal data. It allows simulating correlated and independent samples. Moreover, normality tests and numeric informations are provided.
#'
#'
#'
#' @param n vector size of samples.
#' @param mean vector of means.
#' @param sigma vector of standard deviations or covariance/correlation matrix.
#' @param coefvar an optional vector of coefficients of variation.
#' @param sigmaSup an optional vector of standard deviations if sigma is a correlation matrix.
#' @param dec number of decimals for observations.
#'
#' @usage
#' generator(n , mean = 0, sigma = 1, coefvar = NULL,
#'     sigmaSup = NULL, dec = 2)
#'
#' @details If \code{mean} or \code{sigma} are not specified it's assumed the default values of \code{0} and \code{1}.
#'
#'     If \code{coefvar} (= \code{sigma}/\code{mean}) is specified, function omits \code{sigma} and \code{sigmaSup}. It's assumed that independent samples are desired.
#'
#'     Number of samples are choosen by taken the longest parameter (\code{n}, \code{mean}, \code{sigma}, \code{coefvar}). Therefore, function \code{rep} is used. Pay attention if vectors don't have same length!
#'
#'     If \code{sigma} is a vector, samples are independent. In other case (\code{sigma} is a matrix), samples are dependent (following information meanst be taken into account: if \code{sigma} is a correlation matrix, \code{sigmaSup} is required).
#'
#'
#' @return
#' List containing the following components for independent (with the same length) and dependent samples:
#'
#' \itemize{
#'
#'  \item \code{Samples}: a data frame containing the samples created.
#'
#'
#'  \item Test normality test for the data (\code{shapiro.test()} for n <= 50 and \code{lillie.test()} in other case).
#'
#' }
#'
#' List containing the following components for independent samples with different lengths:
#'
#' \itemize{
#'
#'  \item \code{X_i} sample number i.
#'
#' }
#'
#'
#' @examples
#' generator(4,0,2)
#'
#' sigma <- matrix(c(1,0.8,0.8,1),nrow = 2, byrow = 2)
#' d <- generator(4,mean = c(1,2),sigma, sigmaSup = 1)
#'
#' generator(10,1,coefvar = c(0.3,0.5))
#'
#' generator(c(10,11,10),c(1,2),coefvar = c(0.3,0.5))
#'
#'
#' @export
generator <- function(n, mean = 0, sigma = 1, coefvar = NULL, sigmaSup = NULL, dec = 2){

  if(is.vector(sigma) && any(sigma <= 0))
    stop("Sigma has some negative or null element")

  else if(!is.null(coefvar) && (any(coefvar <= 0) || any(coefvar > 1)))
    stop("Coefficient of variation outbound")


  else if (is.null(coefvar)){
    if(all( length(n) == 1 , length(mean) == 1 ,length(sigma) == 1)){

      x <- stats::rnorm(n, mean, sigma)
      scale <- (x - mean(x)) / stats::sd(x)    # Standardized, to get exactly N(0,1).
      w <- round(scale * sigma + mean, dec)

      if(n <= 50)
        test <- stats::shapiro.test(w) # Shapiro-Wilk's Test.
      else
        test <- nortest::lillie.test(w) # Kolmogorov-Smirnov's Test with Lilliefors correction

      summ <- psych::describe(w)

      l <- list(Sample = w , Test = test, Descriptives = summ)

    }
    else{

      if(is.vector(sigma)){ #Sigma is vector -> independent samples

        k <- max(length(n), length(mean), length(sigma))
        mean <- rep(mean,k)[c(1:k)]
        n <- rep(n,k)[c(1:k)]
        sigma <- rep(sigma,k)[c(1:k)]
        # We force all vectors to have the same length.
        #Therefore, if the user enters vectors with different lengths, it's taken into account the rep function.

        if(length(unique(n)) > 1){ # Different sample sizes are specified.
          l <- list()
          indx <- 1:k

          for(i in indx){
            x <- stats::rnorm(n[i], mean[i] , sigma[i])
            scale <- ( x - mean(x) ) / stats::sd(x)   # Standardized, to get exactly N(0,1).
            w <- round(scale * sigma[i] + mean[i], dec)
            l[[i]] <- w
          }

          namesl <- c() # Empty vector with names for list elements.

          for(i in indx){
            names <- paste("X",i,sep = "") # Name the elements of the list according to the variable number.
            namesl <- c(namesl,names) # Add each variable.
          }

          names(l) <- namesl

        }
        else{ # Vector of sample sizes is the same

          n <- n[1]
          Sigma <- diag(sigma^2)
          x <- data.frame(MASS::mvrnorm(n,mean,Sigma, empirical = TRUE))

          for(i in 1:length(x)){
            x[,i] <- (x[,i]- mean(x[,i])) / stats::sd(x[,i])
            x[,i] <- round(x[,i]*sigma[i] + mean[i], dec)
          } # Observations with specified decimals.

          if (n > 7){
            test <- MVN::mvn(x) # Henze-Zirkler multivariate normality test

            l <- list(Samples = x, Test = test, Descriptives = test$Descriptives)
          }
          else{
            test <- "multivariate normality test available for n>7"
            summ <- psych::describe(x)

            l <- list(Samples = x, Test = test, Descriptives = summ)
          }

        }
      }

      else{ #Sigma is a correlation/covariance matrix-> Dependent samples.

        if(!is.corrmatrix(sigma) && !is.covmatrix(sigma) )
          stop("Please, enter a correlation or covariance matrix")

        else if(is.corrmatrix(sigma) && is.null(sigmaSup))
          stop("Correlation matrix given. Please, enter a support vector of standard deviations (sigmaSup).")

        else if((is.covmatrix(sigma) || is.corrmatrix(sigma)) && length(unique(n)) > 1)
          stop("Sigma is a nonzero covariance matrix. Correlated samples must have same dimension")

        else{

            k1 <- dim(sigma)[1] #Number of rows of Sigma

            k  <- max(length(n), length(mean) , k1)
            mean <- rep(mean,k)[c(1:k)]
            n  <- n[1]

            if(is.corrmatrix(sigma))
              sigma <- mCorrCov(sigma,sigmaSup)

            x  <- data.frame(MASS::mvrnorm(n,mean,sigma,empirical = TRUE))

            for(i in 1:length(x)){
              x[,i] <- (x[,i]- mean(x[,i])) / stats::sd(x[,i])
              x[,i] <- round(x[,i]*sqrt(diag(sigma)[i]) + mean[i], dec)
            } # Observations with specified decimals.

            if (n > 7){
              test <- MVN::mvn(x) # Henze-Zirkler multivariate normality test

              l <- list(Samples = x, Test = test, Descriptives = test$Descriptives)
            }
            else{
              test <- "Multivariate normality test available for n>7"
              summ <- psych::describe(x)

              l <- list(Samples = x, Test = test, Descriptives = summ)
            }
        }
      }
    }
  }

  else{ # coefvar (coefficient of variation) is given.

    if(any(mean == 0)){
      stop("Vector of means cannot be equal to 0") # By definition.
    }

    else {

      if(all(length(n) == 1, length(mean) == 1, length(coefvar) == 1)){ # Univariate distribution.

        sigma <- abs(mean)*coefvar
        x <- stats::rnorm(n, mean, sigma)
        scale <- ( x - mean(x) ) / stats::sd(x)
        w <- round(scale * sigma + mean, dec)

        if(n <= 50)
          test <- stats::shapiro.test(w) #Shapiro-Wilk's Test.
        else
          test <- nortest::lillie.test(w) # Kolmogorov-Smirnov's Test with Lilliefors correction.

        summ <- psych::describe(w)

        l <- list(Sample = w, Test = test, Descriptives = summ)
      }
      else{ # multivariate distribution.

        k <- max(length(n), length(mean), length(coefvar))
        mean <- rep(mean,k)[c(1:k)]
        coefvar <- rep(coefvar,k)[c(1:k)]
        sigma <- abs(mean)*coefvar
        n <- rep(n,k)[c(1:k)]

        if( length(unique(n)) == 1 ){ # Vectors has same size.

          n <- n[1]
          x <- data.frame(MASS::mvrnorm(n,mean,Sigma = diag(sigma^2), empirical = TRUE))

          for(i in 1:length(x)){
            x[,i] <- (x[,i]- mean(x[,i])) / stats::sd(x[,i])
            x[,i] <- round(x[,i]*sigma[i] + mean[i], dec)
          } # Observations with specified decimals.

          if (n > 7){
            test <- MVN::mvn(x) # Henze-Zirkler multivariate normality test

            l <- list(Samples = x, Test = test, Descriptives = test$Descriptives)
          }
          else{
            test <- "Multivariate normality test available for n>7"
            summ <- psych::describe(x)

            l <- list(Samples = x, Test = test, Descriptives = summ)
          }

        }
        else{ # Vectors has different sizes.

            l <- list()
            indx <- 1:k

            for(i in indx){
              x <- stats::rnorm(n[i], mean[i] , sigma[i])
              scale <- ( x - mean(x) ) / stats::sd(x)
              w <- round(scale * sigma[i] + mean[i], dec)
              l[[i]] <- w
            }

            namesl <- c()

            for(i in indx){
              names <- paste("X",i,sep = "")
              namesl <- c(namesl,names)
            }

            names(l) <- namesl

        }
      }
    }
  }
  return(l)
}

# https://es.mathworks.com/matlabcentral/fileexchange/17931-hzmvntest
