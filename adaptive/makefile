#make file for element error analysis

include ../make.inc

adapt: adapt.f90 error.o getnewnodes.o newmeshgen.o tarea.o read.table.o
	$(FC) $(OPS) adapt.f90 -o adapt *.o

error.o: error.f90
	$(FC) $(OPS) -c error.f90

getnewnodes.o: getnewnodes.f90
	$(FC) $(OPS) -c getnewnodes.f90

newmeshgen.o: newmeshgen.f90
	$(FC) $(OPS) -c newmeshgen.f90

tarea.o: tarea.f90
	$(FC) $(OPS) -c tarea.f90

read.table.o: read.table.f90
	$(FC) -c read.table.f90
