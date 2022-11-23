#!/bin/bash

parallel --jobs 12  "Rscript run_estimation_single_mapbiomas_window.R --row {1} --col {2} --grassland_as_forest --combine_other_non_forest --subsample 0.01 --use_md_as_initial_values_for_em --n_random_starts_md 1 &> log_{1}_{2}_use_md_as_initial_values_for_em.txt" ::: {50001..60001..1000} ::: {1..80001..1000}
