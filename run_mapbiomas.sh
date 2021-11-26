#!/bin/bash


parallel --jobs 8  "Rscript explore_mapbiomas.R --row {1} --col {2} --grassland_as_forest --combine_other_non_forest --subsample 0.04 --skip_ml_if_md_is_diag_dominant > log_{1}_{2}.txt" ::: {70000..80000..500} ::: {1000..8000..500} 

