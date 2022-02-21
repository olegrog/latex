SHELL = /bin/bash
PROJECT = $(shell git ls-files | cut -d/ -f1 | grep \.tex$$ | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib sections/* images/*

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)
	rm -f $(PROJECT)Notes.bib

.PHONY: clean zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,spl,bbl}

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
