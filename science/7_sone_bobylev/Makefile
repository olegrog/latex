SHELL = /bin/bash
PROJECT = sone_bobylev
FIGURES = tikz/geometry.pdf
SUBDIRS = contours vs_kn vs_kn2 tikz
SOURCES = $(PROJECT).tex $(PROJECT).bib $(FIGURES)

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

