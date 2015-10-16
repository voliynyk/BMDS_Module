include ../Makefile.conf

all: Hill

Hill: Hill.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o hillfncs.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a
	$(LN) -o hill Hill.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o hillfncs.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist -lassist $(LMINGW)

	$(STRIP) hill$(BINEXT)
	$(CP) hill$(BINEXT) $(BMDS_HOME)/bin/
#	rm hill.exe

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe hill$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.~?~ *~ *.sav

