FROM r-base

RUN apt-get update
RUN yes | apt-get install gdal-bin libgdal-dev

RUN R -e "install.packages(c('data.table', 'ggplot2', 'parallel', 'raster', 'rgdal', 'Rsolnp', 'stargazer', 'optparse'))"
