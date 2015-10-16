# These plots are not currently built or used: 10quanln 10quanqd 15hybrid

DIRS=	continuous 00Hill 00poly 05nctr 05Nlogist 05raivr \
	10gammhit 10logist 10multista 10probit  10weibull 10cancer 10DichoHill

PLOT_EXE= 00expo.exe 00power.exe 00Hill.exe 00poly.exe \
	  05Nlogist.exe 05nctr.exe 05raivr.exe \
	  10weibull.exe 10logist.exe 10probit.exe 10gammhit.exe \
	  10multista.exe  10cancer.exe 10DichoHill.exe

# These targets are not currently used: 15hybrid.exe, 10quanln.exe, 10quanqd.exe

.PHONY: $(PLOT_EXE)
all: $(PLOT_EXE)

10weibull.exe:
	(cd 10weibull;make -f../rMakefile OBJS=weibull_plot.o EXEBASE=10weibull)

10logist.exe:
	(cd 10logist;make -f../rMakefile OBJS=logist_plot.o EXEBASE=10logist)

10probit.exe:
	(cd 10probit;make -f../rMakefile OBJS=probit_plot.o EXEBASE=10probit)

05Nlogist.exe:
	(cd 05Nlogist;make -f../rMakefile OBJS=nlogist_plot.o EXEBASE=05Nlogist)

05nctr.exe:
	(cd 05nctr;make -f../rMakefile OBJS=nctr_plot.o EXEBASE=05nctr)

05raivr.exe:
	(cd 05raivr;make -f../rMakefile OBJS=RaiVR_plot.o EXEBASE=05raivr)

10gammhit.exe:
	(cd 10gammhit;make -f../rMakefile OBJS=gamm_plot.o EXEBASE=10gammhit)

10multista.exe:
	(cd 10multista;make -f../rMakefile OBJS=multista_plot.o EXEBASE=10multista)

00Hill.exe:
	(cd 00Hill;make -f../rMakefile OBJS=00Hill.o EXEBASE=00Hill)

00expo.exe 00power.exe :
	(cd continuous && $(MAKE) $@)

00poly.exe:
	(cd 00poly;make -f../rMakefile OBJS=Poly_plot.o EXEBASE=00poly)

10cancer.exe:
	(cd 10cancer;make -f../rMakefile OBJS=cancer_plot.o EXEBASE=10cancer)

10DichoHill.exe:
	(cd 10DichoHill;make -f../rMakefile OBJS=dichohill_plot.o EXEBASE=10DichoHill)

TAGS clean realclean:
	@for d in $(DIRS); do \
          (cd $$d && $(MAKE) -f ../rMakefile $@) || exit 1; \
	done

# These targets are not currently used
15hybrid.exe:
	(cd 15hybrid;make -f../rMakefile OBJS=Hybrid_plot.o EXEBASE=15hybrid)

10quanln.exe:
	(cd 10quanln;make -f../rMakefile OBJS=quanln_plot.o EXEBASE=10quanln)

10quanqd.exe:
	(cd 10quanqd;make -f../rMakefile OBJS=quanqd_plot.o EXEBASE=10quanqd)
