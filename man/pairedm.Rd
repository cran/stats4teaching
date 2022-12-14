% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairedm.R
\name{pairedm}
\alias{pairedm}
\title{Paired measures (T-Test & Regression)}
\usage{
pairedm(n, mean = 0, sigma = 1, coefvar = NULL,
        rho = NULL, alternative = c("two.sided", "less", "greater"),
        delta = 0, conf.level = 0.95, dec = 2,
        random = FALSE)
}
\arguments{
\item{n}{size of each sample.}

\item{mean}{vector of means.}

\item{sigma}{vector of standard deviations.}

\item{coefvar}{an optional vector of coefficients of variation.}

\item{rho}{Pearson correlation coefficient (optional). If \code{rho} = \code{NULL} a random covariance matrix is generated by \code{genPositiveDefMat()}.}

\item{alternative}{a character string specifying the alternative hypothesis for T-Test. Must be one of ``two.sided`` (default), ``greater`` or ``less``. Can be specified just the initial letter.}

\item{delta}{true value of the difference in means.}

\item{conf.level}{confidence level for interval in T-Test.}

\item{dec}{number of decimals for observations.}

\item{random}{a logical a logical indicating whether you want a random covariance/variance matrix.}
}
\value{
List containing the following components :

\itemize{

 \item \code{Data}: a data frame containing the samples created.


 \item \code{Model}: linear regression model.

 \item \code{T.Test}: a t-test for the samples.

}
}
\description{
Generates two paired measures. It provides T-test and a simple linear regression model for generated data.
}
\details{
If \code{random} = TRUE, \code{rho} is omitted and \code{sigma} is taken as range for variances of the covariance matrix.
}
\examples{

pairedm(10, mean = c(10,2), sigma = c(1.2,0.7), rho = 0.5, alternative = "g")
pairedm(15, mean =c(1,2), coefvar = 0.1, random = TRUE)

}
\seealso{
\code{[clusterGeneration::genpositiveDefMat()]}
}
