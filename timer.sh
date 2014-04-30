#!/bin/bash
echo "" > timer_output.txt
SIZES=(10 100 200 500 1000 1500 2000 4000)
for i in ${SIZES[@]}; do
    ./apsp -n ${i} -csv >> timer_output.txt
done

