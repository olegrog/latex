SHELL = /bin/bash
PROJECT = couette
FIGURES = tikz/geometry.pdf
SUBDIRS = integrated profiles exclusions gamma distrib_f
SOURCES = $(PROJECT).tex $(PROJECT).bib $(FIGURES) appendix.tex

all: $(PROJECT).pdf

tikz/geometry.pdf: tikz/geometry.tex
	cd tikz && latex geometry > /dev/null

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latexmk -pdf $(PROJECT)

.PHONY: $(SUBDIRS) clean

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{bbl}

