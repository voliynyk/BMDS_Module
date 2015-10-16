include ../Makefile.conf

all: logist

logist: logist.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o logist logist.o  -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) logist$(BINEXT)
	$(CP) logist$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav logist$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.~?~ *~ *.sav
