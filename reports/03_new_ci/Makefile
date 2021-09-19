SHELL = /bin/bash
PROJECT = $(patsubst %.tex,%,$(wildcard *.tex))
SOURCES = $(PROJECT).tex
SUBDIRS = $(wildcard */.)

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES) $(SUBDIRS)
	latexmk -pdf $(PROJECT)

.PHONY: clean $(SUBDIRS)

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{spl,bbl}

