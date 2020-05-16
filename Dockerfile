FROM r-base

RUN R -e "install.packages(c('data.table', 'parallel', 'Rsolnp', 'stargazer'))"
