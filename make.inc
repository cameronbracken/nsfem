#this is the make include file for the navier-stokes adaptive mesh simulation

#FC=gfortran
#OPS=-funroll-loops -finline-functions -O3 -fopenmp -ftree-vectorize

#FC=g95
#OPS=-O3

FC=ifort
OPS=-O3 -openmp
