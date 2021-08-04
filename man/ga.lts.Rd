\name{ga.lts}
\alias{ga.lts}
\title{
Function for estimating the LTS (Least Trimmed Squares) regression parameters using genetic algorithms.
}
\description{
This function estimates the LTS (Least Trimmed Squares) regression parameters using genetic algorithms. LTS is a robust regression
estimator with high breakdown property. LTS is a resistant estimator even the number of outliers is up to half of the data. However,
calculating LTS estimator is computaionally expensive. ga.lts() uses evolutionary algorithms (genetic algorithms by default, optionally differential 
evolution) to construct a basic subset and iterates C-steps as defined in Rousseeuw and van-Driessen (2006). Despite lower time efficiency of the ga.lts(),
estimations have lower mean square errors, as a result of lower biases and lower variances. 
}
\usage{
ga.lts(formula, h = NULL, iters = 2, popsize = 50, lower, upper, 
		csteps = 2, method = "ga", verbose = FALSE)
}
\arguments{
  \item{formula}{
Dependent ~ Independents style formula as same in lm() and glm().
}
  \item{h}{
User defined variable to define the majority of the data. Default is floor(n/2)+floor((p+1)/2) where n is the number of observations and p is the number of parameters to estimate
}
  \item{iters}{
Number of generations of the evolutionary algorithm. This variable can be kept larger if the precision is more important than the computation time. Default is 2.
}
  \item{popsize}{
Number of candidates (chromosomes) in the population of evolutionary algorithm. Default is 50.
}
  \item{lower}{
Lower bound for the initial estimates of the parameters.
}
  \item{upper}{
Upper bound for the initial estimates of the parameters.
}
  \item{csteps}{
Number of C-steps to be performed for each candidate solution. Default is 2.
}
  \item{method}{
A string variable for the evolutionary algorithm. 'ga' for the genetic algorithms and 'de' for the differential evolution. Default is 'ga'
}
  \item{verbose}{
A boolean variable for printing the current status of algorithms to screen. Default is FALSE.
}
}

\value{
\item{coefficients}{A vector of estimated parameters}
\item{crit}{LTS criterion value for the reported coefficients}
\item{method}{Name of the method used in the optimization process}
}
\references{
Rousseeuw, P. J., van Driessen, K. (2006). Computing LTS Regression for Large Data Sets ,Data Mining and Knowledge Discovery, 12, 29-45.

Satman, M.,H. (2012). A Genetic Algorithm Based Modification on the LTS Algorithm for Large Data Sets, Communications in Statistics - Simulation and Computation, Vol 41, Issue 5, pp. 644-652.
}
\author{
Mehmet Hakan Satman
}


\examples{
# Data generating process
x1 <- rnorm(100)
x2 <- rnorm(100)
e <- rnorm(100)

# Setting betas to 5
y <- 5 + 5 * x1 + 5 * x2 + e

# Contaminate the data on the dimension of X's randomly
# This is the maximum contamination rate that the LTS can cope with.
outlyings <- sample(1:100, 48)
x1[outlyings] <- 10
x2[outlyings] <- 10

# Estimating LTS with ga (Default optimization method)
lts <- ga.lts(y ~ x1 + x2, popsize = 40, iters = 2, lower = -20, upper = 20)
print(lts)


# Estimating LTS with differential evolution
lts <- ga.lts(y ~ x1 + x2, popsize = 40, iters = 2, lower = -20, upper = 20, method = "de")
print(lts)
}


