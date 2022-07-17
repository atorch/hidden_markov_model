# Hidden Markov Model for Land Cover

This repo includes R code for simulating data from a hidden Markov model for land cover
and estimating parameters using Expectation-Maximization (EM) and Minimum Distance (MD) estimators.

## Citations

Please cite "Improving Estimates of Transitions from Satellite Data: A Hidden Markov Model Approach"
by Adrian L. Torchiana, Ted Rosenbaum, Paul T. Scott, and Eduardo Souza-Rodrigues.

## R Simulations

We run our code in a Docker container in order so that our environment is reproducible. Start by running
```bash
sudo docker build -f ~/hidden_markov_model/Dockerfile ~/hidden_markov_model --tag=hidden_markov_model
sudo docker run -it -v ~/hidden_markov_model:/home/hidden_markov_model hidden_markov_model bash
cd /home/hidden_markov_model/src/
Rscript simulation_simple.R
```

The code in `simulation_simple.R` runs simple simulations that are used to confirm
that the code is working correctly. Its output is not used in the paper.

To reproduce the baseline results in section `D.2 Baseline Results`, run
```bash
Rscript simulation_baseline.R
```

To reproduce section `D.4 Spatial Dependence and Serial Correlation`, run
```bash
Rscript simulation_spatial_corr.R --n_simulations 100
Rscript simulation_spatial_corr.R --n_simulations 100 --z_constant_over_time
```