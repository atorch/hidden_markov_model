# Hidden Markov Model for Land Cover

This repo includes R code for simulating data from a hidden Markov model for land cover
and estimating parameters using Expectation-Maximization (EM) and Minimum Distance (MD) estimators.

## Citations

Please cite "Improving Estimates of Transitions from Satellite Data: A Hidden Markov Model Approach"
by Adrian L. Torchiana, Ted Rosenbaum, Paul T. Scott, and Eduardo Souza-Rodrigues.

## R Simulations

```bash
sudo docker build -f ~/hidden_markov_model/Dockerfile ~/hidden_markov_model --tag=hidden_markov_model
sudo docker run -it -v ~/hidden_markov_model:/home/hidden_markov_model hidden_markov_model bash
cd /home/hidden_markov_model/src/
Rscript simulation_simple.R
Rscript simulation_spatial_corr.R
```