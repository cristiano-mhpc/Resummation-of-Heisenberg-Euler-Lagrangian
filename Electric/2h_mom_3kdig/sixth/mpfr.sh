#!/bin/bash


rm  sixth 

mpicxx -I/usr/local/include/boost/ -o sixth  sixth.cpp  -L/usr/local/lib -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun --oversubscribe -np 1 ./sixth 
