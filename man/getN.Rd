% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learningCurve.R
\name{getN}
\alias{getN}
\title{Function to get the required N given by the accuracy and
the learning curve model}
\usage{
getN(res, acc = 0.9)
}
\arguments{
\item{res}{model results returned by \code{learning_curve} function}

\item{acc}{accuracy that are quired}
}
\value{
sample size that are required
}
\description{
Function to get the required N given by the accuracy and
the learning curve model
}
\examples{
set.seed(2019)
n <- seq(20, 10000, 100)
accMat <- do.call(cbind, lapply(1:length(n), function(i){
tmp_n <- rep(n[i], 50)
y <- -2/(tmp_n^0.8) + 0.95 + rnorm(length(tmp_n), 0, 0.02)
}))
res <- learningCurve(accMat = accMat, n)
N <- getN(res, acc = 0.9)

}
