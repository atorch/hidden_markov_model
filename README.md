# Hidden Markov Model for Land Cover

This repo includes R code for
- simulating land cover data generated by a hidden Markov model, and
- estimating parameters using maximum likelihood (ML using the expectation-maximization algorithm) and minimum distance (MD) estimators.

## Citations

Please cite "Improving Estimates of Transitions from Satellite Data: A Hidden Markov Model Approach"
by Adrian L. Torchiana,
[Ted Rosenbaum](https://www.tedrosenbaum.org/),
[Paul T. Scott](http://ptscott.com/), and
[Eduardo Souza-Rodrigues](https://souza-rodrigues.economics.utoronto.ca/).

## How to Run the Code

We run our code in a Docker container so that our environment is reproducible.
Start by [installing Docker](https://docs.docker.com/engine/install/),
[cloning this repo](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository),
and changing directory into the root of the repo:

```bash
cd hidden_markov_model
```

At this point your working directory should be the root of the repo, and if you `ls` you should see
several .R files as well as the README.md file that you are currently reading.

Start by running
```bash
sudo docker build -f Dockerfile . --tag=hidden_markov_model
```

The `build` command above will use the [Dockerfile](Dockerfile) to create an image with R version 4.2.1 and
[several R packages](install_packages.R) installed.
We can now use that image to run scripts:

```bash
sudo docker run -it -v /home/$USER/hidden_markov_model:/home/hidden_markov_model hidden_markov_model bash
Rscript simulation_simple.R
```

Note that the `-v` in the `docker run` command above is mounting a volume, and assumes that you cloned this repo
into `/home/$USER/hidden_markov_model`. If you cloned into some other location,
replace `/home/$USER/hidden_markov_model` with the correct path.

The code in [simulation_simple.R](simulation_simple.R) runs simple simulations that are used to
confirm that the estimation functions are working correctly. Its output is not used in the paper.

All R scripts should be run inside the container, i.e. following the `docker run` command as above.
Most of the simulations will save output files to the [output](output) directory.

### Baseline Results

To reproduce section `D.2 Baseline Results`, run
```bash
Rscript simulation_baseline.R --n_simulations 100
Rscript describe_simulation_baseline.R --simulation_date yyyy-mm-dd
```
where `yyyy-mm-dd` is the date on which you ran [simulation_baseline.R](simulation_baseline.R).

### Spatial Dependence and Serial Correlation

To reproduce section `D.4 Spatial Dependence and Serial Correlation`, run
```bash
Rscript simulation_spatial_corr.R --n_simulations 100
Rscript simulation_spatial_corr.R --n_simulations 100 --z_constant_over_time
```

### Embrapa Validation

The code for section `6 Validation Exercise Using Land Cover Data` is
[fit_embrapa_validation_hmm.R](fit_embrapa_validation_hmm.R). While we have
made the Embrapa validation code public, the input dataset used in this section
is private and is not included in the repo.

### Atlantic Forest

The code for section `7 Empirical Exercise: Carbon Stocks in the Atlantic Forest`
has several steps:

```bash
mkdir -p atlantic_forest_output
./run_mapbiomas_all_windows.sh
Rscript combine_mapbiomas_estimates.R
./run_viterbi.sh
Rscript atlantic_forest_carbon_stock.R
```

Like section 6, the code for section 7 is public but the input files (the "raw" mapbiomas rasters) are private.