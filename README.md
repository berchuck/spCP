# spCP
Implements a spatially varying change point model with unique intercepts, slopes, variance intercepts and slopes, and change points at each location. Inference is within the Bayesian setting using Markov chain Monte Carlo (MCMC). The response variable can be modeled as Gaussian (no nugget), probit or Tobit link and the five spatially varying parameter are modeled jointly using a multivariate conditional autoregressive (MCAR) prior. The MCAR is a unique process that allows for a dissimilarity metric to dictate the local spatial dependencies. See `vignette('spCP-example')` for usage. Furthermore, the details of the package can be found in the corresponding paper on arXiv by Berchuck et al (2018): "A spatially varying change points model for monitoring glaucoma progression using visual field data", https://arxiv.org/abs/1811.11038v1. 

Finally, we present Video B.1 from the manuscript. At each location on the visual field is the posterior probability of an observed change point presented throughout follow-up and predicted one half-year into the future. By presenting the change-points in video format, the present and future pattern of progression becomes clear.

![Output sample](man/figures/VideoB1.gif)


