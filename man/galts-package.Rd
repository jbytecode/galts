\name{galts-package}
\alias{galts-package}
\alias{galts}
\docType{package}
\title{
Genetic algorithms and C-steps based LTS (Least Trimmed Squares) estimation
}
\description{
This package includes the ga.lts function that estimates LTS (Least Trimmed Squares) parameters using genetic algorithms and C-steps. ga.lts()  constructs a genetic algorithm to form a basic subset and iterates C-steps as defined in Rousseeuw and van-Driessen (2006) to calculate the cost value of the LTS criterion. OLS (Ordinary 
Least Squares) regression is known to be sensitive to outliers. A single outlying observation can change the values of estimated parameters. LTS is a resistant estimator even the number of outliers is up to half of the data. This package is for estimating the LTS parameters with lower bias and variance in a reasonable time. 
}
\details{
\tabular{ll}{
Package: \tab galts\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2011-02-04\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Mehmet Hakan Satman

Maintainer: Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>
}
\references{
Rousseeuw, P. J., van Driessen, K. (2006). Computing LTS Regression for Large Data Sets ,Data Mining and Knowledge Discovery, 12, 29-45.

Satman, M.,H. (2012). A Genetic Algorithm Based Modification on the LTS Algorithm for Large Data Sets, Communications in Statistics - Simulation and Computation, Vol 41, Issue 5, pp. 644-652.
}
\keyword{ package }

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


#Estimating LTS with differential evolution
lts <- ga.lts(y ~ x1 + x2, popsize = 40, iters = 2, lower = -20, upper = 20, method = "de")
print(lts)
}

