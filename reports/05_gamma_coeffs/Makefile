PROJECT = gamma_coeffs
SUBDIRS = pics
SOURCES = $(PROJECT).tex

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

