\name{medmad}
\alias{medmad}
\title{
Function for detecting regression outliers
}
\description{
A method for detecting regression outliers. 
}
\usage{
medmad(formula, h=NULL, csteps=20)
}
\arguments{
  \item{formula}{
Dependent ~ Independents style formula as same in lm() and glm().
}
  \item{h}{
User defined variable to define the majority of the data. Default is floor(n/2)+floor((p+1)/2) where n is the number of observations and p is the number of parameters to estimate
}
  \item{csteps}{
Number of C-steps to be performed for each candidate solution. Default is 2.
}
}

\value{
\item{coefficients}{A vector of estimated parameters}
\item{crit}{LTS criterion value for the reported coefficients}
\item{residuals}{Calculated residuals from the final estimate of model}
}
\author{
Mehmet Hakan Satman
}


\examples{
	n <- 100
	x1 <- rnorm(n,0,10)
	x2 <- rnorm(n,0,10)
	x3 <- rnorm(n,0,10)
	x4 <- rnorm(n,0,10)
	e <- rnorm(n)
	x <- cbind(1, x1, x2, x3, x4)
	p <- 5
	betas <- rep(5,p)
	c <- 0.20
	h <- n - n*c
	y <- 5 + 5*x1 + 5*x2 + 5*x3 + 5*x4 + e
	x1[(h + 1):n] <- rnorm(n-h, 100, 10)
	x2[(h + 1):n] <- rnorm(n-h, 100, 10)
	x3[(h + 1):n] <- rnorm(n-h, 100, 10)
	x4[(h + 1):n] <- rnorm(n-h, 100, 10)

	mm <- medmad(formula = y ~ x1 + x2 + x3 + x4, csteps = 10)

}


