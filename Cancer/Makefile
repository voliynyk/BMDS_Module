include ../Makefile.conf

all: cancer

cancer: cancer.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o cancer cancer.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(STRIP) cancer$(BINEXT)
	$(CP) cancer$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav cancer$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.*~ *~ *.sav


