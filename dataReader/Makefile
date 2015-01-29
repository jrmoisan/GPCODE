netcdflib =libnetcdfreader.a

SRCS =  NetcdfReader_module.f90 MXDReader_module.f90 \
	COLORReader_module.f90 main.f90
	
OBJS =$(SRCS:.f90=.o)	

CC = cc
CFLAGS = -O
FC = /opt/local/bin/mpif90-openmpi-gcc49
F90 = /opt/local/bin/mpif90-openmpi-gcc49

#FC = gfortran
#F90 = gfortran

LDFLAGS = -L/opt/local/lib \
          -I/opt/local/include
LIBS = -lnetcdff
all:: main.x
main.x: $(OBJS)
	$(F90) -o main.x  $(OBJS) $(LDFLAGS) $(LIBS)
clean:
	rm -f $(netcdflib) $(OBJS) *.mod *.o main.x

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) -c $< $(LDFLAGS) $(LIBS)

MXDReader_module.o: NetcdfReader_module.o
COLORReader_module.o: NetcdfReader_module.o
main.o: MXDReader_module.o COLORReader_module.o
