include ../Makefile.conf

all: RaiVR

RaiVR: RaiVR.o initialparms.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o raivr RaiVR.o initialparms.o  -L$(BMDS_HOME)/Assist -lassist $(GSL_LIB) $(LMINGW) 
	$(STRIP) raivr$(BINEXT)
	$(CP) raivr$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe *.sav raivr$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *.~*~ *~ *.sav


