SHELL = /bin/bash
PROJECT = $(shell ls *.tex | cut -f 1 -d '.')
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
	mkdir -p $(PICS)/nonuniform
	cd $(PICS) && ln -sf ../$(COUETTE)/*.pdf .
	cd $(PICS) && ln -sf ../$(COUETTE)/acceleration/*.pdf .
	cd $(PICS)/nonuniform && ln -sf ../../$(COUETTE)/nonuniform/*.pdf .
	cd $(PICS) && ln -sf ../$(POISEUILLE)/*.pdf

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
