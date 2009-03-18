#Top level makefile for adaptive mesh generaton

include make.inc

all: initialmesh sim element_error vectorplot/vp triangle

initialmesh:
	(cd meshgen; make)

sim: 
	(cd nssim; make)
	
element_error: 
	(cd adaptive; make)

vectorplot/vp: vectorplot/vector_plot.f90
	(cd vectorplot; $(FC) vector_plot.f90 -o vp)

triangle: triangle/triangle
	(cd triangle; make)

clean:
	rm datafiles/mesh* datafiles/*.out time.out summary.out \
	triangle/triangle triangle/showme  \
	vectorplot/vp vectorplot/plotinfo.out \
	nssim/subs/*.o nssim/slap/*.o nssim/lib/*.a nssim/nssim \
	adaptive/adapt adaptive/*.o \
	meshgen/meshgen meshgen/*.o meshgen/holefile meshgen/domain

