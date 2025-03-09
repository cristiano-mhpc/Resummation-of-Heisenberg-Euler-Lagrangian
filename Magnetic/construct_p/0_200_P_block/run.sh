#!/bin/bash 
rm -rf build



cmake -S . -B build 
cmake --build build 

mpirun --oversubscribe -np 5  build/auto
