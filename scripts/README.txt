file descriptions for this folder

stan_utility.R: contains function that prints out diagnostic statistics for model fits. Assess effective sample size per iteration (n_eff/iter), RHat, percent divergences, tree depth, and energy Bayesian fraction of missing information (E-BFMI)

write_submodels_D1.Rmd: contains Stan sub-models for Dataset 1 and Dataset 3 (here, \beta_1 is fixed at zero, which is not the case for Dataset 2 models).

write_submodels_D2.Rmd: contains Stan sub-models for Dataset 2. Same models as in write_submodels_D1.R, except \beta_1 is a parameter.

figures.Rmd: contains code for replicating all figures, including SI figures. For the code to run correctly, the full models must be fit; for figure 3, the entire scripts fit_full_models_D[1--3].R must be run first in order to make the three panels.

fit_full_models_D1.Rmd: fits the full model to Dataset 1

fit_full_models_D1.Rmd: same as above but for Dataset 2

fit_full_models_D3.Rmd: same as above but for Dataset 3. Notably, there is a two-step model fitting process, where the posterior from fitting a subset of the data is used as the prior in the second step.

support_functions.R: self explanatory. Mostly useful for the truncated normal distribution (rnormt).

submodel_LOOCV.Rmd: applies Pareto smoothed importance sampling leave one out cross validation (PSIS-LOOCV), a method of model comparisons, to the sub-models of each dataset. The Rmd file produces Tables S2-S4 (as .txt files).

submodel_diagnostics: for each dataset and each submodel, we look at the diagnostic printout (via stan_utility.R) and the posterior contraction. Template code for pair plots is provided. Template code for the shinystan interface is provided.



