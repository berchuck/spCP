# spCP
Implements a spatially varying change point model with unique intercepts, slopes, variance intercepts and slopes, and change points at each location. Inference is within the Bayesian setting using Markov chain Monte Carlo (MCMC). The response variable can be modeled as Gaussian (no nugget), probit or Tobit link and the five spatially varying parameter are modeled jointly using a multivariate conditional autoregressive (MCAR) prior. The MCAR is a unique process that allows for a dissimilarity metric to dictate the local spatial dependencies. Full details of the package can be found in the accompanying vignette.

See `vignette('spCP-example')` for usage.
