include ../Makefile.conf

all: nctr

nctr: Nctr.o initialparms.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o nctr Nctr.o initialparms.o  -L$(BMDS_HOME)/Assist -lassist $(GSL_LIB) $(LMINGW) 

	$(STRIP) nctr$(BINEXT)
	$(CP) nctr$(BINEXT) $(BMDS_HOME)/bin/
#	rm nctr$(BINEXT)


TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe *.sav nctr$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~ *.sav
