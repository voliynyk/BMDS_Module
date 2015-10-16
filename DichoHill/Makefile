include ../Makefile.conf

all: DichoHill

DichoHill: DichoHill.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o DichoHill DichoHill.o  -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) DichoHill$(BINEXT)
	$(CP) DichoHill$(BINEXT) $(BMDS_HOME)/bin/


TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe make.sav DichoHill$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.~?~ *~ *.sav
