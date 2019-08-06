# Hidden Markov Model Simulations

This repo includes both Python and R code for simulating data from a hidden Markov model
and estimating parameters using Expectation-Maximization (EM) and Minimum Distance (MD) estimators.

## Python Simulations

Install docker and then run the following:

```bash
export DOCKER_TAG=hidden_markov_model_docker
sudo docker build ~/hidden_markov_model --tag=$DOCKER_TAG
sudo docker run -it -v ~/hidden_markov_model:/home/hidden_markov_model $DOCKER_TAG bash
cd /home/hidden_markov_model
python src/simulation.py --simulation simple_two_states
python src/simulation.py --simulation simple_three_states
python src/simulation.py --simulation time_varying
```

## R Simulations

TODO Put the .R code in its own directory, document here.

# References

The expectation-maximization (EM) and Viterbi code in [src/estimation.py](src/estimation.py)
is based on https://web.math.princeton.edu/~rvan/orf557/hmm080728.pdf,
but modified to allow for time-varying transition probabilities and panel datasets.