SHELL = /bin/bash
PROJECT = $(shell git ls-tree head --name-only | grep tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT)

.PHONY: clean zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
