#!/bin/bash

echo "Starting first group with grassland as forest"
for colVal in {23000..27000..500}
do
    Rscript explore_mapbiomas.R --row 90000 --col $colVal --grassland_as_forest --subsample 0.05 &
done
wait

echo "Starting first group with separate grassland"
for colVal in {23000..27000..500}
do
    Rscript explore_mapbiomas.R --row 90000 --col $colVal --subsample 0.05 &
done
wait
