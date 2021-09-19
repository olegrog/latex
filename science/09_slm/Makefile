SHELL = /bin/bash
PROJECT = $(shell ls *.tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib
SUBDIRS = single_track experiment experiment/contours
PICS = _pics

all: pics $(PROJECT).pdf

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latexmk -pdf $(PROJECT)

.PHONY: pics clean zip

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}

zip: $(SOURCES)
	zip -r $(PROJECT).zip $(SOURCES) $(PROJECT).bbl $(PICS)
