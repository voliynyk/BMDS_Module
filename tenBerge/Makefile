include ../Makefile.conf

all: tenberge

tenberge: tenberge.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o tenberge tenberge.o -L$(BMDS_HOME)/Assist -lassist $(LMINGW)
	$(STRIP) tenberge$(BINEXT)
	$(CP) tenberge$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	rm -f *.o *.exe make.sav tenberge$(BINEXT)

realclean: clean
	rm -f *.out *.002 *.plt *.~?~ *.sav log.txt
