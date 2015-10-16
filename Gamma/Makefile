include ../Makefile.conf

all: Gamma

Gamma: Gamma.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o gamma Gamma.o  -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) gamma$(BINEXT)
	$(CP) gamma$(BINEXT) $(BMDS_HOME)/bin/


TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav gamma$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.~?~ *~ *.sav
