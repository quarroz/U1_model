#!/bin/bash

method='U1'  # Can be BDF,scipyRK or RK4

for j in `seq 0 3`;
do
    for k in `seq 0 10`;
        do
            qsub /lpt/jquarroz/python/U1_python/script/script$j.sh $k $method
        done
done

