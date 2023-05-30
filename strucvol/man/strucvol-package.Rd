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
  the extended (structural) ARSV model. The multivariate model is implemented using a quasi-maximum likelihood approach (see Harvey, Shephard and Ruiz (1994)). The simulation functions make it possible to analyze the joint behaviour of a stock and a market index using a multivariate structural stochastic volatility model.
  Misspecification tests include: 1. Regression-based LM-type test of a null ARSV model vs. a model with a
  dependent variable in the state equation. 2. Regression-based LM-type test for (remaining) volatility asymmetry ("the leverage effect").
  3. Log likelihood test for testing nested models against each other, presumably ARSV against an extended model. The LM-type tests are conducted by transforming the stochastic volatility model to its log-garch representation (see Asai, 1998) and resolving the resulting distributional misspecification by utilizing the robust tests in Wooldridge (1990). See Back (2023) for details.
}
\author{
Alexander Back

Maintainer: Alexander Back
}
\references{
  Back, Alexander (2023). Structurally motivated stochastic volatility models. Unpublished manuscript.
  
  Sandmann, G. and Koopman, S. J., (1998). Estimation of stochastic volatility models via Monte Carlo maximum likelihood. Journal of Econometrics, 87, issue 2, p. 271-301. https://EconPapers.repec.org/RePEc:eee:econom:v:87:y:1998:i:2:p:271-301.
  
  Asai, M. (1998). A new method to estimate stochastic volatility models: A log-GARCH approach. Journal of the Japan Statistical Society, 28, issue 1, p. 101-114.
  
  Harvey, A., Ruiz, E., & Shephard, N. (1994). Multivariate stochastic variance models. The Review of Economic Studies, 61, issue 2, p. 247-264. https://www.tandfonline.com/doi/abs/10.1198/073500105000000306.
  
  Wooldridge, J. M. (1990). A Unified Approach to Robust, Regression-Based Specification Tests. Econometric Theory, 6, issue 1, p. 17–43. http://www.jstor.org/stable/3532053.
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
