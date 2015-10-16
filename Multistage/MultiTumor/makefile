include ../Makefile.conf

all: ms_combo

ms_combo: ms_combo.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o ms_combo ms_combo.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(STRIP) ms_combo$(BINEXT)
#	$(CP) ms_combo$(BINEXT) $/c:/oakridge/methanol/revised/
#	rm ms_combo$(BINEXT)

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav ms_combo$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.sav


