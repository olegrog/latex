SHELL = /bin/bash
PROJECT = sone_bobylev
FIGURES = tikz/geometry.pdf
SUBDIRS = contours
SOURCES = $(PROJECT).tex $(PROJECT).bib $(FIGURES)

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

tikz/geometry.pdf: tikz/geometry.tex
	cd tikz && latex geometry > log.latex

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

