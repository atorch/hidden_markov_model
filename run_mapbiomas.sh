#!/bin/bash

echo "Starting first group with grassland as forest"

for rowVal in {501..90001..500}
do
    echo "Row Val = $rowVal"
    parallel --jobs 8  "Rscript explore_mapbiomas.R --row $rowVal --col {} --grassland_as_forest --subsample 0.02" ::: {1..80001..500}
done



# echo "Starting first group with separate grassland"
# for colVal in {23000..27000..500}
# do
#     Rscript explore_mapbiomas.R --row 90000 --col $colVal --subsample 0.05 &
# done
# wait
