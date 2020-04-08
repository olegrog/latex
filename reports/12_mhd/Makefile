PROJECT = $(shell git ls-tree head --name-only | grep tex | cut -f 1 -d '.')
SOURCES = $(PROJECT).tex $(PROJECT).bib

all: $(PROJECT).pdf

$(PROJECT).pdf: $(SOURCES)
	latexmk -pdf $(PROJECT) --shell-escape

.PHONY: clean

clean:
	latexmk -C $(PROJECT)
	rm -f $(PROJECT).{run.xml,bbl,gnuploterrors}
	rm -rf gnuplottex
