#make file for initial mesh generator

include ../make.inc

meshgen: regmeshgen.f90 boundary_nodes.o interior_nodes.o obstruction.o
	$(FC) $(OPS) regmeshgen.f90 -o meshgen *.o
	
boundary_nodes.o: boundary_nodes.f90
	$(FC) $(OPS) -c boundary_nodes.f90

interior_nodes.o: interior_nodes.f90
	$(FC) $(OPS) -c interior_nodes.f90
	
obstruction.o: obstruction.f90
	$(FC) $(OPS) -c obstruction.f90
