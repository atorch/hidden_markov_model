# Hidden Markov Model Simulations

Install docker and then run the following:

```bash
export DOCKER_TAG=hidden_markov_model_docker
sudo docker build ~/hidden_markov_model --tag=$DOCKER_TAG
sudo docker run -it -v ~/hidden_markov_model:/home/hidden_markov_model $DOCKER_TAG bash
cd /home/hidden_markov_model
python src/simulation.py
```
