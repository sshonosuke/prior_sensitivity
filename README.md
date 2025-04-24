# Prior sensitivity analysis without model re-fit 

This repository provides R code to ilustrate the Monte Carlo algorithm for prior sensitivity analysis without model re-fit, as proposed by the following paper.

Sugasawa, S. (2024). Prior Sensitivity Analysis without Model Re-fit. [arXiv:2409.19729](https://arxiv.org/abs/2409.19729)

The repository includes the following 6 files.

- `Normal-model.R`: Example of normal distribution
- `BB-model.R`: Example of binomial-beta model
- `GP-model.R`: Example of Gaussian process regression model
- `BB-model-para1.stan` and `BB-model-para2.stan`: Stan codes for binnomial-beta model with different parameterization, used in `BB-model.R`
- `GP-model.stan`: Stan code for Gaussian process regression used in `GP-model.R`



