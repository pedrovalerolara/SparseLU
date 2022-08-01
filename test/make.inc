CC		    = gcc

MKL_LIB     = $(MKLROOT)/lib/intel64  
MKL_INC     = $(MKLROOT)/include

LIBS        = -L$(MKL_LIB) -lm -lmkl_intel_lp64 -lmkl_core -lmkl_sequential 
INC         = -I$(MKL_INC) 

#FLAGS   	= -O3 -Wall --ompss --instrumentation 
FLAGS    	= -DLASs_WITH_MKL -O3 -Wall -fopenmp 
