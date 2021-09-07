SHELL = /bin/bash
PROJECT = couette
FIGURES = tikz/geometry.pdf
SUBDIRS = integrated
SOURCES = $(PROJECT).tex $(PROJECT).bib $(FIGURES) appendix.tex

all: $(PROJECT).pdf

tikz/geometry.pdf: tikz/geometry.tex
	cd tikz && latex geometry > /dev/null

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latex $(PROJECT)
	biber $(PROJECT)
	latex $(PROJECT)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	rm -f $(PROJECT).{aux,toc,log,out,run.xml,bbl,bcf,blg,pdf}

