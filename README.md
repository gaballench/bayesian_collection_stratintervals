# bayesian_collection_stratintervals

Code for the paper on Bayesian stratintervals (Ballen, 2025).

## Structure of this repository

The following directories contain the code for replicating the simulation results in the paper. Each of the following set of directories will be placed into a directory called `lambda_*` where `*` is the specific value of lambda in used, e.g. `lambda_zero` = $\lambda = 0.0$, and `lambda_minushalf` = $\lambda = -0.5$.

- `bayesian_inference_largenumbers`: This directory contains 1000 simulations for testing the Bayesian model.
- `lambda_on_thetas`: The simulations on the effect of different values of $\lambda$ on $\theta_1,\theta_2$
- `single_posteriors`: The posterior sampling from specific posteriors of a fixed, increasing sample size.
- `likelihood_surfaces`: The plots showing the undesirable behaviour of maximum likelihood estimation.

These directories contain the code and data for the empirical examples in the paper.

- `example2_barracudas`: Here we have both code and data for the example on the origination age of the Barracudas.
- `example1_cerrejon`: Example using the Palynomorphs of the Cerrejón formation for showing inference of stratigraphic intervals as well as conflation

# Literature cited

Ballen, G.A. 2025. A flexible Bayesian method for estimating stratigraphic intervals and their co-occurrence in time. BioRxiv. XX:XX-XX
