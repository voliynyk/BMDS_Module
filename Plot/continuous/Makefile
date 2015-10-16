include ../../Makefile.conf
BMDS_HOME = ../..

LN	  = $(CC)
LNFLAGS	  = -L$(BMDS_HOME)/Assist -lassist $(LMINGW) -lm
#INCLUDES  = -I.

all: 00expo.exe 00power.exe

00expo.exe : expoplot.o commonplot.o
	$(LN) -o $@ $^ $(LNFLAGS)
#	$(STRIP) $@
	$(CP) $@ $(BMDS_HOME)/bin/ 

00power.exe : powplot.o commonplot.o
	$(LN) -o $@ $^ $(LNFLAGS)
#	$(STRIP) $@
	$(CP) $@ $(BMDS_HOME)/bin/ 

commonplot.o : commonplot.c
	$(CC) -c $(CFLAGS) $<

#contplot.dll : commonplot.o
#	$(LN) $(DLL_LN_FLAGS) -o $@ $^

# Implicit rules
.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

# Miscellaneous targets

TAGS:
	etags *.c --include=$(BMDS_HOME)/Assist/TAGS

clean:
	rm -f *.o *.exe *.dll *.lib *.a

realclean: clean
	rm -f *~ \#*\#
