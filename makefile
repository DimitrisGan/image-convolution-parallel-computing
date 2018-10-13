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

mpi:
	$(CMPI) $(CFLAGS) -o $(MPDIR)mpi_exe $(mpi_source) -lm

#-----------------------------------------------------------------------------------

open_mp:
	$(CMPI) $(CFLAGS) $(CMPIFLAGS) -o $(OPENMPDIR)open_mp_exe $(open_mp_source) -lm

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
