include ../Makefile.conf

all: multistage cancer

multistage: multistage.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o multistage multistage.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(STRIP) multistage$(BINEXT)
	$(CP) multistage$(BINEXT) $(BMDS_HOME)/bin/

cancer: cancer.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o cancer cancer.o LoadCommBloc.o getmle.o getcl.o donlp2usrfc.o -L$(BMDS_HOME)/Assist  -lassist $(LMINGW) 
	$(STRIP) cancer$(BINEXT)
	$(CP) cancer$(BINEXT) $(BMDS_HOME)/bin/

cancer.o: multistage.c
	$(CC) -o cancer.o -c multistage.c -DBMDS_CANCER_MODEL $(CFLAGS)

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav multistage$(BINEXT) cancer$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.sav


