#!/bin/bash

parallel --jobs 4  "Rscript run_viterbi.R --row {1} --col {2} --subsample 0.01 --width_in_pixels 1000 --raster_year 2017 &> viterbi_log_{1}_{2}_2017.txt" ::: {1..100000..1000} ::: {1..80001..1000}

