\name{strucvol-package}
\alias{strucvol-package}
\alias{strucvol}
\docType{package}
\title{
  Estimate and test structural stochastic volatility models.
}
\description{
  This package implements routines for estimation and misspecification testing 
  in the context of structural stochastic volatility models. In addition, it 
  can be used to model other stochastic volatility models with dependent variables
  in the state equation of the underlying state space system.
}
\details{
  The Monte Carlo maximum likelihood
  method in Sandmann and Koopmann (1996) is implemented to estimate the standard ARSV model and 
  the extended model. The multivariate model is implemented using a quasi-maximum likelihood approach.
  Misspecification tests include: 1. Regression-based LM-type test of a null ARSV model vs. a model with a
  dependent variable in the state equation. 2. Regression-based LM-type test for (remaining) volatility asymmetry ("the leverage effect").
  3. Log likelihood test for testing nested models against each other, presumably ARSV against an extended model.
}
\author{
Alexander Back

Maintainer: Alexander Back
}
\references{
  Back (2023),
  Sandmann and Koopmann (1996),
  Wooldridge (1988).
}
\keyword{Stochastic volatility, Monte Carlo likelihood, state space models, misspecification testing}
\seealso{
...
}
\examples{
  \dontrun{
     ## Try using the function fitsv() on the first column in the dataframe "df".
  }
}
