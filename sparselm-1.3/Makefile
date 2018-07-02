#
# GNU Makefile for the sparseLM library
#

CC=gcc
GCCFLAGS=-O3 -Wall -ffast-math -funroll-all-loops #-g #-fomit-frame-pointer # -msse3 # -fprofile-arcs -ftest-coverage # -pg # GCC-specific optimization flags
MACH_CFLAGS=-march=pentiumpro -malign-double # machine specific gcc flags, in this case for the pentium pro
#MACH_CFLAGS=-march=native -mtune=native

SRCS=splm.c splm_crsm.c splm_ccsm.c splm_hessian.c splm_time.c splm_stm.c splm_misc.c
OBJS=splm.o splm_crsm.o splm_ccsm.o splm_hessian.o splm_time.o splm_stm.o splm_misc.o

AR=ar
RANLIB=ranlib
MAKE=make


############ Start of configuration options ############

## read in available sparse solvers
include solvers.mk

## header file paths
SUITESPARSE_PATH=/home/lourakis/SuiteSparse-3.4.0

CSPARSE_INC=-I$(SUITESPARSE_PATH)/CSparse/Include
LDL_INC=-I$(SUITESPARSE_PATH)/LDL/Include
CHOLMOD_INC=-I$(SUITESPARSE_PATH)/CHOLMOD/Include
PARDISO_INC=
UMFPACK_INC=-I$(SUITESPARSE_PATH)/UMFPACK/Include
SuperLU_INC=-I/usr/local/include/slu
TAUCS_INC=-I/usr/local/include/
SPOOLES_INC=-I/home/lourakis/spooles.2.2/
MUMPS_INC=-I/home/lourakis/MUMPS_4.8.3/include/ -I/home/lourakis/MUMPS_4.8.3/libseq
MKL_INC=-I/opt/intel/mkl/8.0.1/include/
UFconfig_INC=-I$(SUITESPARSE_PATH)/UFconfig/
AMD_INC=-I$(SUITESPARSE_PATH)/AMD/Include/
METIS_INC=

######## No changes should be necessary beyond this point ########


ifeq ($(strip $(HAVE_CSPARSE)), yes)
  CFL:=$(CFL) -DHAVE_CSPARSE
  INCLUDES:=$(INCLUDES) $(CSPARSE_INC)
  OBJS:=$(OBJS) splm_csparse.o
  SRCS:=$(SRCS) splm_csparse.c
  NEED_UFconfig=no
endif

ifeq ($(strip $(HAVE_LDL)), yes)
  CFL:=$(CFL) -DHAVE_LDL
  INCLUDES:=$(INCLUDES) $(LDL_INC)
  OBJS:=$(OBJS) splm_ldlp.o
  SRCS:=$(SRCS) splm_ldlp.c
  NEED_UFconfig=yes
  NEED_AMD=yes
endif

ifeq ($(strip $(HAVE_CHOLMOD)), yes)
  CFL:=$(CFL) -DHAVE_CHOLMOD
  INCLUDES:=$(INCLUDES) $(CHOLMOD_INC)
  OBJS:=$(OBJS) splm_cholmod.o
  SRCS:=$(SRCS) splm_cholmod.c
  NEED_UFconfig=yes
endif

ifeq ($(strip $(HAVE_MA57)), yes)
  CFL:=$(CFL) -DHAVE_MA57
  OBJS:=$(OBJS) splm_ma57.o ma57.o
  SRCS:=$(SRCS) splm_ma57.c ma57.c
  NEED_METIS=yes
endif

ifeq ($(strip $(HAVE_MA47)), yes)
  CFL:=$(CFL) -DHAVE_MA47
  OBJS:=$(OBJS) splm_ma47.o ma47.o
  SRCS:=$(SRCS) splm_ma47.c ma47.c
endif

ifeq ($(strip $(HAVE_MA27)), yes)
  CFL:=$(CFL) -DHAVE_MA27
  OBJS:=$(OBJS) splm_ma27.o ma27.o
  SRCS:=$(SRCS) splm_ma27.c ma27.c
endif

ifeq ($(strip $(HAVE_PARDISO)), yes)
  CFL:=$(CFL) -DHAVE_PARDISO
  INCLUDES:=$(INCLUDES) $(PARDISO_INC)
  OBJS:=$(OBJS) splm_pardiso.o
  SRCS:=$(SRCS) splm_pardiso.c
endif

ifeq ($(strip $(HAVE_SuperLU)), yes)
  CFL:=$(CFL) -DHAVE_SuperLU
  INCLUDES:=$(INCLUDES) $(SuperLU_INC)
  OBJS:=$(OBJS) splm_superlu.o
  SRCS:=$(SRCS) splm_superlu.c
endif

ifeq ($(strip $(HAVE_UMFPACK)), yes)
  CFL:=$(CFL) -DHAVE_UMFPACK
  INCLUDES:=$(INCLUDES) $(UMFPACK_INC)
  OBJS:=$(OBJS) splm_umfpack.o
  SRCS:=$(SRCS) splm_umfpack.c
  NEED_UFconfig=yes
endif

ifeq ($(strip $(HAVE_TAUCS)), yes)
  CFL:=$(CFL) -DHAVE_TAUCS
  INCLUDES:=$(INCLUDES) $(TAUCS_INC)
  OBJS:=$(OBJS) splm_taucs.o
  SRCS:=$(SRCS) splm_taucs.c
  NEED_METIS=yes
endif

ifeq ($(strip $(HAVE_SPOOLES)), yes)
  CFL:=$(CFL) -DHAVE_SPOOLES
  INCLUDES:=$(INCLUDES) $(SPOOLES_INC)
  OBJS:=$(OBJS) splm_spooles.o
  SRCS:=$(SRCS) splm_spooles.c
endif

ifeq ($(strip $(HAVE_MUMPS)), yes)
  CFL:=$(CFL) -DHAVE_MUMPS
  INCLUDES:=$(INCLUDES) $(MUMPS_INC)
  OBJS:=$(OBJS) splm_mumps.o
  SRCS:=$(SRCS) splm_mumps.c
endif

ifeq ($(strip $(HAVE_MA77)), yes)
  CFL:=$(CFL) -DHAVE_MA77
  OBJS:=$(OBJS) splm_ma77.o
  SRCS:=$(SRCS) splm_ma77.c
  NEED_METIS=yes
endif

ifeq ($(strip $(HAVE_MKL)), yes)
  CFL:=$(CFL) -DHAVE_MKL
  INCLUDES:=$(INCLUDES) $(MKL_INC)
  OBJS:=$(OBJS) splm_dss.o
  SRCS:=$(SRCS) splm_dss.c
endif


ifeq ($(strip $(NEED_UFconfig)), yes)
  #CFL:=$(CFL)
  INCLUDES:=$(INCLUDES) $(UFconfig_INC)
endif

######## ordering options ########
ifeq ($(strip $(NEED_AMD)), yes)
  INCLUDES:=$(INCLUDES) $(AMD_INC)
endif

ifeq ($(strip $(NEED_METIS)), yes)
  CFL:=$(CFL) -DUSE_METIS
  INCLUDES:=$(INCLUDES) $(METIS_INC)
endif


CFLAGS=$(INCLUDES) $(GCCFLAGS) $(MACH_CFLAGS) $(CFL) -DUSE_ATA_QUADS
# DSS: -I/opt/intel/mkl/8.0.1/include/

## dependencies

all: libsplm.a dem

libsplm.a: $(OBJS)
	$(AR) crv libsplm.a $(OBJS)
	$(RANLIB) libsplm.a

splm.o: splm.h splm_priv.h compiler.h
splm_crsm.o: splm.h compiler.h
splm_ccsm.o: splm.h compiler.h
splm_hessian.o: splm.h splm_priv.h compiler.h
splm_time.o: splm.h compiler.h
splm_stm.o: splm.h
splm_misc.o: splm.h splm_priv.h
splm_ldlp.o: splm.h splm_priv.h compiler.h
splm_ma47.o: splm.h splm_priv.h compiler.h
splm_ma27.o: splm.h splm_priv.h compiler.h
splm_pardiso.o: splm.h splm_priv.h
splm_superlu.o: splm.h splm_priv.h
splm_taucs.o: splm.h splm_priv.h
splm_csparse.o: compiler.h splm_priv.h
splm_ma57.o: splm.h compiler.h splm_priv.h
splm_cholmod.o: splm.h compiler.h splm_priv.h
splm_umfpack.o: splm.h compiler.h splm_priv.h
splm_spooles.o: splm.h compiler.h splm_priv.h
splm_mumps.o: splm.h compiler.h splm_priv.h
splm_dss.o: splm.h splm_priv.h
splm_ma77.o: splm.h compiler.h splm_priv.h

dem:
	$(MAKE) -f Makefile_demo

clean:
	@rm -f $(OBJS)
	$(MAKE) -f Makefile_demo clean

cleanall: clean
	@rm -f libsplm.a
	$(MAKE) -f Makefile_demo cleanall

depend:
	makedepend -f cond.make $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

