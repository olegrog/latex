SHELL = /bin/bash
PROJECT = $(shell git ls-tree head --name-only | grep tex$$ | cut -f 1 -d '.')
SUBDIRS = validation
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latexmk -pdf $(PROJECT) --shell-escape

.PHONY: clean zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl,gnuploterrors}
	rm -rf gnuplottex

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
