% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tauCens.R
\name{tauCens}
\alias{tauCens}
\title{Estimator of Kendall's tau for censored data}
\usage{
tauCens(formula, data = NULL, method = "adjusted")
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~
operator, and the terms on the right.  The response must be a
survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.}

\item{data}{a data.frame containing the variables in the model.}

\item{method}{which estimator to use. The non-parametric estimator from Hougaard (2000)
or the naive one where non-fully observed pairs are left out.}
}
\value{
non-parametric estimate of Kendall's tau and parametric estimates from three different
frailty models.
}
\description{
Estimator of Kendall's tau for censored data
}
\details{
Kendall's tau is a rank based measure of dependence. Note that this estimator is biased towards zero (and so is the estimator for the variance).
}
\references{
Hougaard, Philip. Analysis of multivariate survival data. Springer Science & Business Media, 2012.
}
\seealso{
tauPar taucpp
}
\author{
Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
}
