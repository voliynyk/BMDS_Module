include ../Makefile.conf

all: Weibull

Weibull: Weibull.o $(BMDS_HOME)/Assist/libassist.a ProfLik.o Makefile
	$(LN) -o weibull Weibull.o ProfLik.o  -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) weibull$(BINEXT)
	$(CP) weibull$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe weibull$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav

