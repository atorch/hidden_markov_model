# Hidden Markov Model Simulations

This repo includes both Python and R code for simulating data from a hidden Markov model
and estimating parameters using Expectation-Maximization (EM) and Minimum Distance (MD) estimators.

## R Simulations

```bash
export DOCKER_TAG_R=hidden_markov_model_r
sudo docker build -f ~/hidden_markov_model/Dockerfile-R ~/hidden_markov_model --tag=$DOCKER_TAG_R
sudo docker run -it -v ~/hidden_markov_model:/home/hidden_markov_model $DOCKER_TAG_R bash
cd /home/hidden_markov_model/src/
Rscript simulation_counties.R
```

# References

The expectation-maximization (EM) and Viterbi code in [src/python/estimation.py](src/python/estimation.py)
is based on https://web.math.princeton.edu/~rvan/orf557/hmm080728.pdf,
with modifications to allow for time-varying transition probabilities and panel datasets.
The most relevant sections in Ramon van Handel's HMM notes are

* Algorithm 3.1: Forward-Backward Algorithm
* Algorithm 3.2: Baum-Welch Algorithm
* Algorithm 3.4: Viterbi Algorithm
* Algorithm 6.1: Concrete EM Algorithm
