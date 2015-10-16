include ../Makefile.conf

CFLAGS := $(CFLAGS) $(GSL_INCLUDE)

all: nlogist

nlogist: Nlogist.o $(BMDS_HOME)/Assist/libassist.a Makefile
	$(LN) -o nlogist Nlogist.o  -L$(BMDS_HOME)/Assist -lassist $(GSL_LIB) $(LMINGW) 
	$(STRIP) nlogist$(BINEXT)
	$(CP) nlogist$(BINEXT) $(BMDS_HOME)/bin/

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	$(RM) *.o *.exe *.sav nlogist$(BINEXT)

realclean: clean
	$(RM) *.out *.002 *.plt *~ *.*~*.sav


