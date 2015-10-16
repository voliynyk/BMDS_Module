include ../Makefile.conf

all: Poly

Poly: Poly.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o polyfncs.o getmle.o getcl.o  donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a 
	$(LN) -o poly Poly.o LoadCommBloc.o A3fcns.o ContFcns.o GETMLEA3.o polyfncs.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) poly$(BINEXT)
	$(CP) poly$(BINEXT) $(BMDS_HOME)/bin/


TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav poly$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav

