#!/bin/bash

for j in `seq 0 3`;
do
    for k in `seq 0 10`;
        do
            python3 /home/jeremie/Thesis/Code/U1_model/U1/U1.py $j $k
        done
done

