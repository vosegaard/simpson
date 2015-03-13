SRC = B0inhom.c OCroutines.c allocation.c auxmath.c averaging.c blockdiag.c \
  cm.c complx.c cryst.c crystdat.c fft.c fidcalc.c ftools.c ham.c iodata.c isotopes.c lbfgs.c \
  main.c matrix.c pthread_barrier_mac.c pulse.c readsys.c relax.c rfprof.c rfshapes.c \
  sim.c simpson.c spinach.c spinsys.c tclcode.c tclutil.c wigner.c
OBJ = $(SRC:.c=.o)

# WINDOWS:
#INCLUDES = -IC:/Tcl/include -I../CBLAS/src -I/usr/local/include -I../fftw3
#LIBRARIES = -lm libfftw3-3.dll libnfft3-0.dll tcl85.dll blas.dll cblas.dll lapack.dll -ldl -lpthread
# Linux (Ubuntu 14.04):
INCLUDES = -I/usr/include/tcl8.5
LIBRARIES = -ltcl8.5 -llapack -lblas -lfftw3 /usr/local/lib/libnfft3.a -lpthread -lm
FLAGS = -DNO_CONST -O3

CC = gcc
RM = rm
TAR = tar

simpson: $(OBJ)
	$(CC) $(FLAGS) $(OBJ) $(LIBRARIES) -o simpson 
.c.o:
	$(CC) $(FLAGS) $(INCLUDES) -c $<
clean:
	$(RM) -f *.o simpson
dist:
	$(TAR) cvzf simpson.tgz *.c *.h simpson.xcodeproj Makefile
