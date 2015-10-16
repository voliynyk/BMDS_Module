include ../Makefile.conf

all: ms_combo

ms_combo: ms_combo.o loadcommbloc.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o ms_combo ms_combo.o loadcommbloc.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(STRIP) ms_combo$(BINEXT)
	$(CP) ms_combo$(BINEXT) $(BMDS_HOME)/bin/
#	rm ms_combo$(BINEXT)

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav ms_combo$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.sav

testcombo: testcombo.f donlp2usrfc.o Makefile ../Makefile.conf
	$(FC) -c testcombo.f $(FFLAGS)
	$(LN) -o testcombo testcombo.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(CP) testcombo$(BINEXT) $(BMDS_HOME)/bin/
