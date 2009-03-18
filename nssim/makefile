#makefile for navier stokes simulation: nssim

include ../make.inc

MKLPATH=/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/32
MKLINCLUDE=/Library/Frameworks/Intel_MKL.framework/Versions/10.0.1.014/include
ACCELERATE=-Wl,-framework -Wl,Accelerate
LIBS=$(MKLPATH)/libmkl_intel.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a $(MKLPATH)/libguide.a -lpthread
#-I$(MKLINCLUDE) -L$(MKLPATH) $(LIBS)

nssim: lib/subs.a lib/libdlap.a main_ns.f90
	$(FC) $(OPS) main_ns.f90 -o nssim lib/subs.a lib/libdlap.a $(ACCELERATE)

####Below here for for navier-stokes simulation subroutines

lib/subs.a: subs/boundary_type.o subs/boundary_conditions.o subs/assemble_jb.o\
            subs/read.table.o subs/jbsubs.o subs/residual_fem.o \
            subs/boundary_adjust_residual.o subs/jacobian_fem.o \
            subs/residual_fem_jb.o subs/jacobian_fem_jb.o
	ar -r lib/subs.a subs/*.o

subs/boundary_type.o: subs/boundary_type.f90
	(cd subs; $(FC) $(OPT) -c boundary_type.f90)

subs/boundary_conditions.o: subs/boundary_conditions.f90
	(cd subs; $(FC) $(OPT) -c boundary_conditions.f90)

subs/assemble_jb.o: subs/assemble_jb.f90
	(cd subs; $(FC) $(OPT) -c assemble_jb.f90)

subs/read.table.o: subs/read.table.f90
	(cd subs; $(FC) $(OPT) -c read.table.f90)

subs/jbsubs.o: subs/jbsubs.f90
	(cd subs; $(FC) $(OPT) -c jbsubs.f90)

subs/residual_fem.o: subs/residual_fem.f90
	(cd subs; $(FC) $(OPT) -c residual_fem.f90)

subs/boundary_adjust_residual.o: subs/boundary_adjust_residual.f90
	(cd subs; $(FC) $(OPT) -c boundary_adjust_residual.f90)

subs/jacobian_fem.o: subs/jacobian_fem.f90
	(cd subs; $(FC) $(OPT) -c jacobian_fem.f90)

subs/residual_fem_jb.o: subs/residual_fem_jb.f90
	(cd subs; $(FC) $(OPT) -c residual_fem_jb.f90)

subs/jacobian_fem_jb.o: subs/jacobian_fem_jb.f90
	(cd subs; $(FC) $(OPT) -c jacobian_fem_jb.f90)

#####Below here for compiling slap library libdlap and interface program slapit

lib/libdlap.a: slap/slapit.o slap/dlap.o
	ar -r lib/libdlap.a slap/*.o

slap/slapit.o: slap/slapit.f90
	(cd slap; $(FC) $(OPT) -c slapit.f90)

slap/dlap.o: slap/dlap.f90
	(cd slap; $(FC) $(OPT) -c dlap.f90)
