#!/bin/bash


parallel --jobs 8  "Rscript explore_mapbiomas.R --row {1} --col {2} --grassland_as_forest --combine_other_non_forest --subsample 0.04 --skip_ml_if_md_is_diag_dominant --n_random_starts_md 1 > log_{1}_{2}.txt" ::: {1..50001..1000} ::: {1..80000..1000} 

