ICC = icc
MPICPP = mpic++
# MPIRUN = mpirun
UPCC = upcc
# UPCRUN = upcrun

# number of processors for MPI code
ifndef NPES
	NPES = 4
endif

# network API for UPC
ifndef NETWK
	NETWK = smp
endif

# mpi flags (specifying number of processors)
MPI_FLAGS = -DNPES=$(NPES)
# upcc flags
UPCC_FLAGS = -network=$(NETWK) -T=$(NPES) -DNPES=$(NPES)


targets = laplace_serial laplace_mpi laplace_upc

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

laplace_serial: laplace_serial.c
	$(ICC) -o $@ laplace_serial.c
laplace_mpi: laplace_mpi.c
	$(MPICPP) $(MPI_FLAGS) -o $@ laplace_mpi.c
laplace_upc: laplace_upc.c 
	$(UPCC) $(UPCC_FLAGS) -o $@ laplace_upc.c 


.PHONY : clean
clean:
	rm -f $(targets) *.stdout *.txt


# DIRS = cpi doublereduce laplace mg
# TOP = $(shell pwd)
# THREADS=2
# icc:
# 	$(ICC) laplace

# clean:
# 	for d in $(DIRS); do\
# 	(cd $$d && $(MAKE) clean)\
# 	done

# runall:
# 	for d in $(DIRS); do\
# 	(cd $$d && $(MAKE) run THREADS=$(THREADS))\
# 	done

