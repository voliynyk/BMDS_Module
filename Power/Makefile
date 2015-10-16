include ../Makefile.conf

all: Power

Power: Power.o ProfLik.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o powfncs.o getmle.o getcl.o getprofile.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a
	$(LN) -o power Power.o ProfLik.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o powfncs.o getmle.o getcl.o getprofile.o donlp2usrfc.o -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) power$(BINEXT)
	$(CP) power$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe power$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav

