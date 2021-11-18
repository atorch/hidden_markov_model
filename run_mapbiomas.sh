#!/bin/bash

echo "Starting first group with grassland as forest"
for colVal in {23000..27000..500}
do
    Rscript explore_mapbiomas.R --row 90000 --col $colVal --grassland_as_forest&
done
wait

echo "Starting first group with separate grassland"
for colVal in {23000..27000..500}
do
    Rscript explore_mapbiomas.R --row 90000 --col $colVal&
done
wait


# echo "Starting second group"
# Rscript explore_mapbiomas.R --row 90200 --col 23000 &
# Rscript explore_mapbiomas.R --row 90200 --col 23200 &
# Rscript explore_mapbiomas.R --row 90200 --col 23400 &
# Rscript explore_mapbiomas.R --row 90200 --col 23600 &
# Rscript explore_mapbiomas.R --row 90200 --col 23800 &
# Rscript explore_mapbiomas.R --row 90200 --col 24000 &
# Rscript explore_mapbiomas.R --row 90200 --col 24200 &
# Rscript explore_mapbiomas.R --row 90200 --col 24400 &
# Rscript explore_mapbiomas.R --row 90200 --col 24600 &
# Rscript explore_mapbiomas.R --row 90200 --col 24800 &
# wait

# echo "Starting third group"
# Rscript explore_mapbiomas.R --row 90400 --col 23000 &
# Rscript explore_mapbiomas.R --row 90400 --col 23200 &
# Rscript explore_mapbiomas.R --row 90400 --col 23400 &
# Rscript explore_mapbiomas.R --row 90400 --col 23600 &
# Rscript explore_mapbiomas.R --row 90400 --col 23800 &
# Rscript explore_mapbiomas.R --row 90400 --col 24000 &
# Rscript explore_mapbiomas.R --row 90400 --col 24200 &
# Rscript explore_mapbiomas.R --row 90400 --col 24400 &
# Rscript explore_mapbiomas.R --row 90400 --col 24600 &
# Rscript explore_mapbiomas.R --row 90400 --col 24800 &
# wait

# echo "Starting fourth group"
# Rscript explore_mapbiomas.R --row 90600 --col 23000 &
# Rscript explore_mapbiomas.R --row 90600 --col 23200 &
# Rscript explore_mapbiomas.R --row 90600 --col 23400 &
# Rscript explore_mapbiomas.R --row 90600 --col 23600 &
# Rscript explore_mapbiomas.R --row 90600 --col 23800 &
# Rscript explore_mapbiomas.R --row 90600 --col 24000 &
# Rscript explore_mapbiomas.R --row 90600 --col 24200 &
# Rscript explore_mapbiomas.R --row 90600 --col 24400 &
# Rscript explore_mapbiomas.R --row 90600 --col 24600 &
# # Rscript explore_mapbiomas.R --row 90600 --col 24800 &
# wait

# echo "Starting fifth group (move to rows north of 90000)"
# Rscript explore_mapbiomas.R --row 89800 --col 23000 &
# Rscript explore_mapbiomas.R --row 89800 --col 23200 &
# Rscript explore_mapbiomas.R --row 89800 --col 23400 &
# Rscript explore_mapbiomas.R --row 89800 --col 23600 &
# Rscript explore_mapbiomas.R --row 89800 --col 23800 &
# Rscript explore_mapbiomas.R --row 89800 --col 24000 &
# Rscript explore_mapbiomas.R --row 89800 --col 24200 &
# Rscript explore_mapbiomas.R --row 89800 --col 24400 &
# Rscript explore_mapbiomas.R --row 89800 --col 24600 &
# Rscript explore_mapbiomas.R --row 89800 --col 24800 &
# wait

# echo "Starting sixth group"
# Rscript explore_mapbiomas.R --row 89600 --col 23000 &
# Rscript explore_mapbiomas.R --row 89600 --col 23200 &
# Rscript explore_mapbiomas.R --row 89600 --col 23400 &
# Rscript explore_mapbiomas.R --row 89600 --col 23600 &
# Rscript explore_mapbiomas.R --row 89600 --col 23800 &
# Rscript explore_mapbiomas.R --row 89600 --col 24000 &
# Rscript explore_mapbiomas.R --row 89600 --col 24200 &
# Rscript explore_mapbiomas.R --row 89600 --col 24400 &
# Rscript explore_mapbiomas.R --row 89600 --col 24600 &
# Rscript explore_mapbiomas.R --row 89600 --col 24800 &
# wait

# echo "Starting seventh group"
# Rscript explore_mapbiomas.R --row 89400 --col 23000 &
# Rscript explore_mapbiomas.R --row 89400 --col 23200 &
# Rscript explore_mapbiomas.R --row 89400 --col 23400 &
# Rscript explore_mapbiomas.R --row 89400 --col 23600 &
# Rscript explore_mapbiomas.R --row 89400 --col 23800 &
# Rscript explore_mapbiomas.R --row 89400 --col 24000 &
# Rscript explore_mapbiomas.R --row 89400 --col 24200 &
# Rscript explore_mapbiomas.R --row 89400 --col 24400 &
# Rscript explore_mapbiomas.R --row 89400 --col 24600 &
# Rscript explore_mapbiomas.R --row 89400 --col 24800 &
# wait

# echo "Starting eight group"
# Rscript explore_mapbiomas.R --row 89200 --col 23000 &
# Rscript explore_mapbiomas.R --row 89200 --col 23200 &
# Rscript explore_mapbiomas.R --row 89200 --col 23400 &
# Rscript explore_mapbiomas.R --row 89200 --col 23600 &
# Rscript explore_mapbiomas.R --row 89200 --col 23800 &
# Rscript explore_mapbiomas.R --row 89200 --col 24000 &
# Rscript explore_mapbiomas.R --row 89200 --col 24200 &
# Rscript explore_mapbiomas.R --row 89200 --col 24400 &
# Rscript explore_mapbiomas.R --row 89200 --col 24600 &
# Rscript explore_mapbiomas.R --row 89200 --col 24800 &
# wait

# echo "Starting ninth group"
# Rscript explore_mapbiomas.R --row 89000 --col 23000 &
# Rscript explore_mapbiomas.R --row 89000 --col 23200 &
# Rscript explore_mapbiomas.R --row 89000 --col 23400 &
# Rscript explore_mapbiomas.R --row 89000 --col 23600 &
# Rscript explore_mapbiomas.R --row 89000 --col 23800 &
# Rscript explore_mapbiomas.R --row 89000 --col 24000 &
# Rscript explore_mapbiomas.R --row 89000 --col 24200 &
# Rscript explore_mapbiomas.R --row 89000 --col 24400 &
# Rscript explore_mapbiomas.R --row 89000 --col 24600 &
# Rscript explore_mapbiomas.R --row 89000 --col 24800 &
# wait

# echo "Starting tenth group"
# Rscript explore_mapbiomas.R --row 88800 --col 23000 &
# Rscript explore_mapbiomas.R --row 88800 --col 23200 &
# Rscript explore_mapbiomas.R --row 88800 --col 23400 &
# Rscript explore_mapbiomas.R --row 88800 --col 23600 &
# Rscript explore_mapbiomas.R --row 88800 --col 23800 &
# Rscript explore_mapbiomas.R --row 88800 --col 24000 &
# Rscript explore_mapbiomas.R --row 88800 --col 24200 &
# Rscript explore_mapbiomas.R --row 88800 --col 24400 &
# Rscript explore_mapbiomas.R --row 88800 --col 24600 &
# Rscript explore_mapbiomas.R --row 88800 --col 24800 &
# wait

# echo "Starting eleventh group"
# Rscript explore_mapbiomas.R --row 88600 --col 23000 &
# Rscript explore_mapbiomas.R --row 88600 --col 23200 &
# Rscript explore_mapbiomas.R --row 88600 --col 23400 &
# Rscript explore_mapbiomas.R --row 88600 --col 23600 &
# Rscript explore_mapbiomas.R --row 88600 --col 23800 &
# Rscript explore_mapbiomas.R --row 88600 --col 24000 &
# Rscript explore_mapbiomas.R --row 88600 --col 24200 &
# Rscript explore_mapbiomas.R --row 88600 --col 24400 &
# Rscript explore_mapbiomas.R --row 88600 --col 24600 &
# Rscript explore_mapbiomas.R --row 88600 --col 24800 &
# wait
