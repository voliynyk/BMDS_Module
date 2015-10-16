# Makefile for Assist
include ../Makefile.conf

OBJ = ERRORPRT.o gettimeofday.o slogf.o in_outfun.o CLOSE_FILES.o\
      allo_memo.o Dmngb.o run_dmngb.o blas/d1mach.o blas/dgemm.o blas/dger.o \
      blas/dscal.o blas/dswap.o blas/dtrsm.o blas/idamax.o blas/xerbla.o \
      lapack/dgesv.o lapack/dgetf2.o lapack/dgetrf.o lapack/dgetrs.o \
      lapack/dlaswp.o lapack/lsame.o lapack/ieeeck.o lapack/ilaenv.o \
      lapack/dgetri.o lapack/dtrti2.o lapack/dtrtri.o lapack/xerbla.o \
      go/fmin.o zeroin.o CIs.o \
      Nested_CI.o N_Goodness.o N_Bootstrap.o SRoI.o \
      donlp2/user_eval.o donlp2/donlp2.o donlp2/testgrad.o  \
      donlp2/solchk.o donlp2/sysdepen.o \
      Goodness.o \
      INITIALIZE_CVG3.o COMPUTE_DTMS3LOGLKH.o COMPUTE_DTMS3FLOGLKH.o \
      COMPUTE_DTMS3CVG_Gc.o COMPUTE_DTMS3CVG_Pc.o COMPUTEDF.o DTMS3ANOVA.o \
      VARSUMCOMP.o DTMS3ANOVAC.o Take_Out_Bounded_Parms.o Get_Linear_Trend.o \
      SWAP_IVECTOR.o FILL_SPECVECTOR.o COUNT_SPECVECTOR.o GETUNKNOWNPARMS.o \
      GETKNOWNPARMS.o FILL_LDMATRIXBYUD.o INITIALIZE_DVECTOR.o INITIALIZE_DMATRIX.o \
      INITIALIZE_UDMATRIX.o MATMPYV.o MATMPYM.o MATMPYV2.o TRANSPOSE.o \
      MATINVS.o INVMAT.o MATMPYM2.o \
      GAMMP.o XGAMMAI_A.o XGAMMAI.o CNORM.o NORM.o RNORM.o CHISQ.o QCHISQ.o qstudt.o \
      PROBABILITY_W.o PROBABILITY_INRANGE.o Sort_2_By_Dose.o SortByLs.o \
      Sort_4_By_Dose.o Get_Names.o path_name.o path_name2.o Slog.o D_Slog.o \
      dcdflib/ipmpar.o dcdflib/dcdflib.o \
      show_version.o C_Goodness.o Extras.o compute_continuous_liks.o 



CFLAGS = -Wall -g -O3 $(MORE_CFLAGS) -DASSIST
## These flags are for g77
##FFLAGS = -g -O3 -fno-automatic -finit-local-zeros -I.
## These flags are for gfortran
FFLAGS = -g -O3 -fno-automatic -I. 

all: libassist.a

libassist.a: $(OBJ) Makefile
	ar rv libassist.a $(OBJ)
	ranlib libassist.a

TAGS:
	etags *.f *.c *.h blas/*.f dcdflib/*.c dcdflib/*.h lapack/*.f

clean:
	$(RM) *.o *.a
	(cd blas; $(RM) *.o)
	(cd lapack; $(RM) *.o)
	(cd dcdflib; $(RM) *.o)
	(cd go; $(RM) *.o)

realclean: clean
	$(RM) *~ *.*~ *out 
	(cd blas; $(RM) *~ *.*~ *out)
	(cd lapack; $(RM) *.~*~ *~ *out)
	(cd dcdflib; $(RM) *.~*~ *~ *out)
	(cd go; $(RM) *.~*~ *~ *out)
