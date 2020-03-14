\name{NPrior-package}
\alias{NPrior-package}
\alias{NPrior}
\docType{package}
\title{
\packageTitle{NPrior}
}
\description{
This package provides a set of functions that implements a general framework of Bayesian shrinkage priors called "Neuronized Priors". 
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{NPrior}
\packageIndices{NPrior}
}
\examples{
library(MASS)
data(Boston)
str(Boston)
attach(Boston)
X = cbind(crim,indus,nox,rm,age,dis,tax,ptratio,black,lstat)
X = scale(X)
y = medv
y = y-mean(y)
n = nrow(X)
p = ncol(X)
#fit = NPrior_run(X, y)
}
