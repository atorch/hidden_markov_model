FROM r-base:4.2.1

RUN apt-get update
RUN yes | apt-get install gdal-bin libgdal-dev

# First we install the remotes package
RUN R -e "install.packages('remotes')"

# Next, we use remotes::install_version to install specific versions of several R packages
COPY install_packages.R .
RUN Rscript install_packages.R

WORKDIR /home/hidden_markov_model
