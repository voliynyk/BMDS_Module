include Makefile.conf

SUBDIRS = Assist Cancer DichoHill Exponential Gamma Hill Logistic ms_combo Multistage NCTR NLogistic Plot Poly Power Probit RaiVR Weibull tenBerge

all TAGS:: 
	@for d in $(SUBDIRS); do\
		(cd $${d} && $(MAKE) $@) || exit 1;\
	done

clean realclean::
#	$(RM) *~
	(cd bin; $(RM) *)
	@for d in $(SUBDIRS); do\
		(cd $${d} && $(MAKE) $@) || exit 1;\
	done

zip : 
	(DIR=$$(basename "$$(pwd)");cd bin; zip ../../$${DIR}_$$(date +%Y%m%d%H%M)_exe.zip *$(BINEXT))

dist:
	make realclean
	(DIR=$$(basename "$$(pwd)");cd ..;zip -r $${DIR}_$$(date +%Y%m%d%H%M)_srs.zip $${DIR})

check:
	(cd TestData && $(MAKE) $@) || exit 1;

