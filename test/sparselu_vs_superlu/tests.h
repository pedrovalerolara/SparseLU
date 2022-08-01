#include "/home/5pv/SummerIntern21/Cameron-SparseLU/SuperLU_MT_3.1/SRC/slu_mt_ddefs.h"
#include "../../include/lass.h"
#include "unistd.h"
#include "stdio.h"
#include "errno.h"
#include "string.h"
#include "fcntl.h"
#include "sys/wait.h"
#include "sys/sysinfo.h"
#include "stdint.h"

extern void test_superlu(int nprocs,
                         struct timeval *start, struct timeval *end,
                         double *flops);

extern void test_slass(int nprocs, char *filename,
                       struct timeval *start, struct timeval *end,
                       double *flops);
