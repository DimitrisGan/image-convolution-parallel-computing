CMPI=mpicc
CFLAGS=-Wall -g
CMPIFLAGS= -fopenmp

MPDIR=mpi/
mpi_source=$(addprefix $(MPDIR),mpi.c)
mpi_outputs=$(addprefix $(MPDIR),output_images/*)

OPENMPDIR=open_mp/
open_mp_source=$(addprefix $(OPENMPDIR),open_mp.c)
open_mp_outputs=$(addprefix $(OPENMPDIR),output_images/*)


.PHONY: all mpi open_mp

all: mpi open_mp

#-------------------------------------------------------------------------

mpi: mpi.o
	$(CMPI) $(CFLAGS) mpi.o -o $(MPDIR)mpi_exe -lm
	rm -f mpi.o

mpi.o: $(mpi_source)
	$(CMPI) $(CFLAGS) -c $(mpi_source)

#-----------------------------------------------------------------------------------

open_mp: open_mp.o
	$(CMPI) $(CFLAGS) $(CMPIFLAGS) open_mp.o -o $(OPENMPDIR)open_mp_exe -lm
	rm -f open_mp.o

open_mp.o: $(open_mp_source)
	$(CMPI) $(CFLAGS) $(CMPIFLAGS) -c $(open_mp_source)

#-----------------------------------------------------------------------------------

.PHONY: clean_all clean_exes clean_mpi clean_open_mp clean_outs clean_mpi_outs clean_open_mp_outs

clean_all: clean_exes clean_outs

#-----------------------------------------------------------------------------------

clean_exes: clean_mpi clean_open_mp

clean_mpi:
	rm -f $(MPDIR)mpi_exe

clean_open_mp:
	rm -f $(OPENMPDIR)open_mp_exe

#-----------------------------------------------------------------------------------

clean_outs: clean_mpi_outs clean_open_mp_outs

clean_mpi_outs:
	rm -r $(mpi_outputs)

clean_open_mp_outs:
	rm -r $(open_mp_outputs)
