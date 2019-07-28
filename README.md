## Requirements
Compilers: icc, mpi++, upcc


## Compile
By default, the parallel versions (laplace_mpi.c and laplace_upc.c) will be compiled for 4-core machine. One can change this configuration by setting "NPES" argument when calling "make":

make NPES=<number of processors>

The compiled executables must be run with the exact number of processors it was compiled for. For example, "make NPES=8" will create executables laplace_mpi and laplace_upc which must later be run with exactly 8 processors.

To compile the UPC version, one needs to supply a code for network API. By default, this is set to "smp" so that the compiled executable can be run on one single node (e.g. personal laptop). One can change the network API by setting "NETWK" argument when calling "make".


## Run the solver
When running the compiled executables (laplace_serial, laplace_mpi, laplace_upc), one can pass the following arguments:
    -q                      if specified, the solver will run in quiet mode (no console output will be printed)
    -m <max iteration>      the max number of iterations the solver can run
    -s <summary file name>  if specified, final outputs will be appended to the file


## Optional
The autorun.sh bash script is written to compile and run the three versions with different settings for processor number and max iterations to compare their performance.