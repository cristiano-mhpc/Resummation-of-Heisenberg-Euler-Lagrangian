#!/bin/bash 
rm -rf build

cmake -S . -B build 
cmake --build build 
inputs=("20" "25" "30" "35" "40" "50" "60" "100" "200" "500")

#Loop through each input and run the binary

for value in "${inputs[@]}"; do
    echo "Running with $value coefficients"
    build/delta "$((value-1))"
    echo "------------"
done
