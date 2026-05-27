# Bayesian Dynamic Models and MCMC Methods

* In this repository, I implement some Markov Chain Monte Carlo (MCMC) algorithms and Bayesian dynamic models using R. The purpose is to study and better understand these methods. The implementations are not expected to be fully optimized, but may serve as a starting point.

## State Space Dynamic Models

### Poisson Local level model
- [poisson_local_level_sim.R](poisson_local_level/poisson_local_level_sim.R): Data simulation.
- [poisson_sin_level_sim.R](poisson_local_level/poisson_sin_level_sim.R): Data simulation with sin level.
- [poisson_local_level_mh_gibbs.R](poisson_local_level/poisson_local_level_mh_gibbs.R): Metropolis within Gibbs (Geweke and Tanikazi, 2001). 
- [poisson_local_level_cwmh.R](poisson_local_level/poisson_local_level_cwmh.R): Component Wise Metropolis Hastings (CWMH) (Montoril, Correia and Migon 2022).
- [poisson_local_level_stan.R](poisson_local_level/poisson_local_level_cwmh.R): Simulation with [Stan](https://mc-stan.org/rstan/)
- [poisson_pol2.stan](poisson_local_level.stan): Stan model

### Poisson Second Order Polynomial Model
- [poisson_pol2_sim.R](poisson_pol2/poisson_local_level_sim.R): Data simulation
- [poisson_pol2_mh_gibbs.R](poisson_pol2/poisson_pol2_mh_gibbs.R): Metropolis within Gibbs (Geweke and Tanikazi, 2001).
- [poisson_pol2_stan.R](poisson_pol2/poisson_local_level_cwmh.R): Simulation with [Stan](https://mc-stan.org/rstan/) (NUTS)
- [poisson_pol2.stan](poisson_pol2/poisson_pol2.stan): Stan model
- [poisson_schnatter.R](poisson_pol2/poisson_schnatter.R): Gibbs with data augmentation (Frühwirth-Schnatter and Wagner, 2006)
 
## Basic MCMC methods
- [rejection_sampling.R](basic_methods/rejection_sampling.R): Simple rejection sampling. Samples from a target distribution (normalization is not required) using the uniform distribution as the envelope function.
- [adaptive_rejection_sampling.R](basic_methods/adaptive_rejection_sampling.R): Adaptive rejection sampling (tangent method) of Gilks and Wild (1992). Requires log-concave target density.
- [metropolis.R](basic_methods/metropolis.R): Random walking metropolis method.
- [augmented_data.R](basic_methods/augmented_data.R): Implementation of the Linkage Model example using Data Augmentation / Substitution algorithm of Tanner and Wong (1987).

## References
- Frühwirth-Schnatter, S., & Wagner, H. (2006). Auxiliary Mixture Sampling for Parameter-Driven Models of Time Series of Counts with Applications to State Space Modelling. Biometrika, 93(4), 827–841.
- Geweke, J., & Tanizaki, H. (2001). Bayesian estimation of state-space models using the Metropolis–Hastings algorithm within Gibbs sampling. Computational statistics & data analysis, 37(2), 151-170
- Gilks, W. R., e P. Wild. “Adaptive Rejection Sampling for Gibbs Sampling”. Applied Statistics 41, n. 2 (1992): 337. https://doi.org/10.2307/2347565.
- Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains and their applications.
- Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of state calculations by fast computing machines. The journal of chemical physics, 21(6), 1087-1092.4
- Montoril, M. H., Correia, L. T., & Migon, H. S. (2021). Bayesian estimation of dynamic weights in Gaussian mixture models. arXiv preprint arXiv:2104.03395.
- Tanner, M. A., & Wong, W. H. (1987). The calculation of posterior distributions by data augmentation. Journal of the American statistical Association, 82(398), 528-540.