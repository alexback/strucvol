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
  method in Sandmann and Koopmann (1998) is implemented to estimate the standard ARSV model and 
  the extended model. The multivariate model is implemented using a quasi-maximum likelihood approach. The simulation functions make it possible to analyze the joint behaviour of a stock and a market index using a multivariate structural stochastic volatility model.
  Misspecification tests include: 1. Regression-based LM-type test of a null ARSV model vs. a model with a
  dependent variable in the state equation. 2. Regression-based LM-type test for (remaining) volatility asymmetry ("the leverage effect").
  3. Log likelihood test for testing nested models against each other, presumably ARSV against an extended model.
}
\author{
Alexander Back

Maintainer: Alexander Back
}
\references{
  Back, Alexander (2023). Structurally motivated stochastic volatility models. Unpublished manuscript.
  
  Sandmann, Gleb and Koopman, Siem Jan, (1998). Estimation of stochastic volatility models via Monte Carlo maximum likelihood. Journal of Econometrics, 87, issue 2, p. 271-301. https://EconPapers.repec.org/RePEc:eee:econom:v:87:y:1998:i:2:p:271-301.
  
  Wooldridge, Jeffrey M. (1990). A Unified Approach to Robust, Regression-Based Specification Tests. Econometric Theory, 6, issue 1, p. 17–43. http://www.jstor.org/stable/3532053.
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
