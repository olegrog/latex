SHELL = /bin/bash
PROJECT = $(shell git ls-tree head --name-only | grep tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib
COUETTE = ../../reports/9_bgk_lattice/couette
POISEUILLE = poiseuille
PICS = _pics

all: pics $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

.PHONY: pics clean zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}

pics:
	mkdir -p $(PICS)
	cd $(PICS) && ln -sf ../$(COUETTE)/*.pdf .
	cd $(PICS) && ln -sf ../$(COUETTE)/acceleration/*.pdf .
	cd $(PICS) && ln -sf ../$(POISEUILLE)/*.pdf

zip: $(SOURCES) highlights.txt
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
