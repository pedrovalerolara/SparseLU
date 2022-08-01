############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:        March 20, 2013  Version 2.1
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _OPENMP

#
#  The name of the libraries to be created/linked to
#
TMGLIB       = libtmglib.a
SUPERLULIB   = /home/5pv/SummerIntern21/Cameron-SparseLU/SuperLU_MT_3.1/lib/libsuperlu_mt_OPENMP.a

BLASDEF   = -DUSE_VENDOR_BLAS
BLASLIB = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
MATHLIB   = -lm
MPLIB     =

BIN = /home/5pv/SummerIntern21/Cameron-SparseLU/workplan-summer-2021-internship-copy/LASs-sparse-npgetrf-and-dense-tile-pivoting/LASs-sparse-npgetrf/test/BIN

#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

# Definitions used by CPP:
PREDEFS         = -D__OPENMP
#

CC           = gcc
CFLAGS       = -gdwarf-2 -gstrict-dwarf -fopenmp -DLASs_WITH_MKL -Wall -g -O3 -m64 -lm -I"${MKLROOT}/include" -I$(S)/include #-lumfpack -lamd -lcholmod -lsuitesparseconfig -lm
#
# There are a few files that should NOT be optimized, which appear in the
# end of SRC/Makefile, taking the flag $(NOOPTS)
# NOOPTS should be defined to explicitly turn off any optimization.
NOOPTS   = -O0

############################################################################
FORTRAN      = gfortran -fopenmp
FFLAGS       =
LOADER       = gcc -fopenmp
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DNoChange
