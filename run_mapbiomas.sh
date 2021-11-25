#!/bin/bash

for rowVal in {80000..91000..500}
do
    for colVal in {16000..26500..500}
    do
	logFile=log_${rowVal}_${colVal}.txt
	echo "Running row $rowVal col $colVal, writing logs to $logFile"
	Rscript explore_mapbiomas.R --row $rowVal --col $colVal --grassland_as_forest --combine_other_non_forest --subsample 0.04 --skip_ml_if_md_is_diag_dominant &> $logFile &
    done
    echo "Sleeping before running next group of jobs"
    sleep 1200
done
