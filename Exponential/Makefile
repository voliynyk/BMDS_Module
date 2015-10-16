include ../Makefile.conf

all: exponential

exponential: exponential.o ProfLik.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o expfncs.o getmle.o getcl.o getprofile.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a
	$(LN) -o exponential exponential.o ProfLik.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o expfncs.o getmle.o getcl.o getprofile.o donlp2usrfc.o -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) exponential$(BINEXT)
	$(CP) exponential$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe exponential$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav

