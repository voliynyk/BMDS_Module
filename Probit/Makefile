include ../Makefile.conf

all: probit

probit: probit.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o probit probit.o  -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) probit$(BINEXT)
	$(CP) probit$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav probit$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav

