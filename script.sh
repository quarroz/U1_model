#!/bin/bash

for j in `seq 0 100`;
do
    qsub /lpt/jquarroz/python/U1/script_error.sh $j
done
