#!/bin/bash

for rowVal in {89000..91000..500}
do
    for colVal in {20000..25500..500}
    do
	logFile=log_${rowVal}_${colVal}.txt
	echo "Running row $rowVal col $colVal, writing logs to $logFile"
	Rscript explore_mapbiomas.R --row $rowVal --col $colVal --grassland_as_forest --combine_other_non_forest --subsample 0.04 &> $logFile &
    done
    echo "Sleeping before running next group of jobs"
    sleep 3600
done
