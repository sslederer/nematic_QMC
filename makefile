#makefile for Linux using the intel ifort compiler
RM = rm -f

OBJS = dqmc.o        greens_function.o  lapack_interface.o  mc_config.o  module_global.o\
       random_gen.o  savemat.o          tests.o             measure.o    boson_meas.o

LIBDIR= "$(shell pwd)"

FLAGS =  -I$(LIBDIR) -inline-forceinline -no-inline-max-size

LDFLAGS = -L/opt/intel/composer_xe_2013/mkl/lib/intel64/ -L/opt/intel/composer_xe_2013/compiler/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -liomp5

F90=ifort -c $(FLAGS)

dqmc: $(OBJS)
	ifort $(FLAGS) $(OBJS) $(LDFLAGS) -o  dqmc

lib: $(OBJS)
	ar cr libdqmc.a $(OBJS)

module_global.o: module_global.f90
	$(F90) $<

random_gen.o: random_gen.f
	$(F90) $<

%.o: %.f90 random_gen.o module_global.o
	$(F90) $<

clean:
	$(RM) *.o dqmc

all: dqmc lib
