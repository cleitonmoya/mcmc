# Bayesian Dynamic Models and MCMC Methods

* In this repository, I implement some Markov Chain Monte Carlo (MCMC) algorithms and Bayesian Dynamic Models using R. The purpose is to study and better understand them. Implementations are not expected to be optimized, but they can serve as a starting point.

- [rejection_sampling.R](rejection_sampling.R): Simple rejection sampling. Samples from a target distribution (normalization is not required) using the uniform distribution as the envelope function.
- [adaptive_rejection_sampling.R](adaptive_rejection_sampling.R): Adaptive rejection sampling (tangent method) of Gilks and Wild (1992). Requires log-concave target density.
- [metropolis.R](metropolis.R): Random walking metropolis method.
- [augmented_data.R](augmented_data.R): Implementation of the Linkage Model example using Data Augmentation / Substitution algorithm of Tanner and Wong (1987).

# State Space Dynamic Models
- [poisson_local_level_naive_gibbs.R](poisson_local_level_naive_gibbs.R): Poisson local level model using naive Gibbs Sampling. 
- [poisson_local_level_cwmh.R](poisson_local_level_cwmh.R): Poisson local level model using the Component Wise Metropolis Hastings (CWMH) proposal of Montoril, Correia and Migon (2022).
- [poisson_local_level_sim.R](poisson_local_level_sim.R): Simulation of a Poisson local level model time series.

## References
- Gilks, W. R., e P. Wild. “Adaptive Rejection Sampling for Gibbs Sampling”. Applied Statistics 41, n. 2 (1992): 337. https://doi.org/10.2307/2347565.
- Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains and their applications.
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of state calculations by fast computing machines. The journal of chemical physics, 21(6), 1087-1092.4
- Montoril, M. H., Correia, L. T., & Migon, H. S. (2021). Bayesian estimation of dynamic weights in Gaussian mixture models. arXiv preprint arXiv:2104.03395.
- Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior distributions by data augmentation. Journal of the American statistical Association, 82(398), 528-540.