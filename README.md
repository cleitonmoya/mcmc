# Markov Chain Monte Carlo (MCMC) Methods

* In this repository, I implement some MCMC algorithms using R. The purpose is to study and better understand them. Implementations are not expected to be optimized, but they can serve as a starting point.

- [rejection_sampling.R]: Simple rejection sampling. Samples from a target distribution (normalization is not required) using the uniform distribution as the envelope function.
- [adaptive_rejection_sampling.R]: Adaptive rejection sampling (tangent method) of Gilks and Wild (1992). Requires log-concave target density.

## References
- Gilks, W. R., e P. Wild. “Adaptive Rejection Sampling for Gibbs Sampling”. Applied Statistics 41, n. 2 (1992): 337. https://doi.org/10.2307/2347565.
