SHELL = /bin/bash
PROJECT = sone_bobylev
SUBDIRS = contours vs_kn vs_kn2 tikz
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

