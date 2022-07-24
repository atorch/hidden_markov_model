FROM r-base:4.2.1

RUN apt-get update
RUN yes | apt-get install gdal-bin libgdal-dev

RUN R -e "install.packages(c('data.table', 'ggplot2', 'parallel', 'raster', 'rgdal', 'rgeos', 'Rsolnp', 'stargazer', 'optparse', 'ngspatial', 'plot.matrix', 'stringr'))"

RUN yes | apt-get install parallel

WORKDIR /home/hidden_markov_model
