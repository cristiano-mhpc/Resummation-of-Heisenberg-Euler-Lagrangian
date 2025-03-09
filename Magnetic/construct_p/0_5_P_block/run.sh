#!/bin/bash 
rm -rf build

cmake -S . -B build 
cmake --build build 

mpirun --oversubscribe -np 1  build/auto
