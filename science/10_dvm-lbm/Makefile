SHELL = /bin/bash
PROJECT = $(shell ls *.tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib
FIGS = ../../../reports/9_bgk_lattice/couette
PICS = pics

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}

pics:
	mkdir -p $(PICS)
	cd $(PICS) && ln -s $(FIGS)/*.pdf .
	cd $(PICS) &&ln -s $(FIGS)/acceleration/*.pdf .

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PICS)/
